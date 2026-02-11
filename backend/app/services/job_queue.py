from __future__ import annotations

import os
import signal
import subprocess
import traceback
from collections import deque
from datetime import datetime, timezone
from threading import Condition, Event, Thread
from typing import Any

from ..core.config import settings
from ..storage.state_manager import JsonStateStore


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _terminate_process_tree(pid: int) -> bool:
    """Best-effort process termination for a running training subprocess."""
    try:
        if os.name == "nt":
            result = subprocess.run(
                ["taskkill", "/PID", str(pid), "/T", "/F"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            return result.returncode == 0
        os.kill(pid, signal.SIGTERM)
        return True
    except (OSError, subprocess.SubprocessError):
        return False


def _normalize_scheduler_mode(value: str | None) -> str:
    mode = (value or "").strip().lower()
    if mode not in {"serial", "parallel"}:
        return "serial"
    return mode


class JobQueue:
    def __init__(self, store: JsonStateStore, runner: Any, auth_service: Any) -> None:
        self.store = store
        self.runner = runner
        self.auth_service = auth_service
        self._stop_event = Event()
        self._workers: list[Thread] = []
        self._cv = Condition()
        self._pending: deque[str] = deque()
        self._running_per_user: dict[str, int] = {}
        self._worker_count = max(1, int(os.getenv("LABFLOW_JOB_WORKERS", "8")))

    def start(self) -> None:
        if self._workers and any(worker.is_alive() for worker in self._workers):
            return
        self._stop_event.clear()
        self._workers = []
        for i in range(self._worker_count):
            worker = Thread(target=self._work_loop, daemon=True, name=f"job-worker-{i}")
            worker.start()
            self._workers.append(worker)

    def stop(self) -> None:
        self._stop_event.set()
        with self._cv:
            self._cv.notify_all()
        for worker in self._workers:
            worker.join(timeout=3)

    def enqueue(self, job_id: str) -> None:
        with self._cv:
            self._pending.append(job_id)
            self._cv.notify_all()

    def cancel(self, job_id: str) -> dict[str, Any] | None:
        job = self.store.get_job(job_id)
        if job is None:
            return None

        status = str(job.get("status") or "")
        if status in {"success", "failed", "canceled"}:
            return job

        if status == "queued":
            with self._cv:
                self._pending = deque(item for item in self._pending if item != job_id)
            return self.store.update_job(
                job_id,
                {
                    "status": "canceled",
                    "ended_at": _now(),
                    "error_msg": "Canceled by user.",
                    "cancel_requested": True,
                },
            )

        patch: dict[str, Any] = {"cancel_requested": True}
        pid = job.get("train_pid")
        if isinstance(pid, int):
            patch["cancel_signal_sent"] = _terminate_process_tree(pid)
        return self.store.update_job(job_id, patch)

    def _owner_key(self, job: dict[str, Any] | None) -> str:
        if not job:
            return "user:anonymous"
        owner = job.get("owner_user_id")
        if isinstance(owner, int):
            return f"user:{owner}"
        if isinstance(owner, str) and owner.strip():
            return f"user:{owner.strip()}"
        return "user:anonymous"

    def _running_limit_for_owner(self, job: dict[str, Any]) -> int:
        owner = job.get("owner_user_id")
        if not isinstance(owner, int):
            return 1
        try:
            mode = _normalize_scheduler_mode(self.auth_service.get_scheduler_mode(owner))
        except Exception:  # noqa: BLE001
            mode = "serial"
        return 3 if mode == "parallel" else 1

    def _reserve_next_job(self) -> tuple[str, dict[str, Any], str] | None:
        with self._cv:
            while not self._stop_event.is_set():
                stale_ids: set[str] = set()
                selected_job_id: str | None = None
                selected_job: dict[str, Any] | None = None
                selected_owner_key = ""

                for candidate_id in self._pending:
                    candidate_job = self.store.get_job(candidate_id)
                    if candidate_job is None:
                        stale_ids.add(candidate_id)
                        continue
                    if candidate_job.get("status") == "canceled":
                        stale_ids.add(candidate_id)
                        continue
                    if candidate_job.get("status") != "queued":
                        stale_ids.add(candidate_id)
                        continue

                    owner_key = self._owner_key(candidate_job)
                    running = self._running_per_user.get(owner_key, 0)
                    limit = self._running_limit_for_owner(candidate_job)
                    if running < limit:
                        selected_job_id = candidate_id
                        selected_job = candidate_job
                        selected_owner_key = owner_key
                        break

                if stale_ids:
                    self._pending = deque(
                        item for item in self._pending if item not in stale_ids
                    )

                if selected_job_id and selected_job:
                    try:
                        self._pending.remove(selected_job_id)
                    except ValueError:
                        continue
                    self._running_per_user[selected_owner_key] = (
                        self._running_per_user.get(selected_owner_key, 0) + 1
                    )
                    return selected_job_id, selected_job, selected_owner_key

                self._cv.wait(timeout=1.0)

        return None

    def _release_owner_slot(self, owner_key: str) -> None:
        with self._cv:
            current = self._running_per_user.get(owner_key, 0)
            if current <= 1:
                self._running_per_user.pop(owner_key, None)
            else:
                self._running_per_user[owner_key] = current - 1
            self._cv.notify_all()

    def _work_loop(self) -> None:
        while not self._stop_event.is_set():
            reserved = self._reserve_next_job()
            if reserved is None:
                continue

            job_id, job, owner_key = reserved
            try:
                if job.get("status") == "canceled":
                    continue

                self.store.update_job(
                    job_id,
                    {"status": "running", "started_at": _now(), "error_msg": None},
                )
                if job["type"] == "train":
                    self._execute_train(job_id, job)
                elif job["type"] == "predict":
                    self._execute_predict(job_id, job)
                else:
                    raise RuntimeError(f"Unsupported job type: {job['type']}")
                latest = self.store.get_job(job_id) or {}
                if latest.get("status") != "canceled":
                    self.store.update_job(job_id, {"status": "success", "ended_at": _now()})
            except Exception as exc:  # noqa: BLE001
                latest = self.store.get_job(job_id) or {}
                if latest.get("cancel_requested"):
                    self.store.update_job(
                        job_id,
                        {
                            "status": "canceled",
                            "ended_at": _now(),
                            "error_msg": "Canceled by user.",
                        },
                    )
                else:
                    self.store.update_job(
                        job_id,
                        {
                            "status": "failed",
                            "ended_at": _now(),
                            "error_msg": f"{exc}\n{traceback.format_exc()}",
                        },
                    )
            finally:
                self._release_owner_slot(owner_key)

    def _execute_train(self, job_id: str, job: dict[str, Any]) -> None:
        dataset = self.store.get_dataset(job["dataset_id"])
        if dataset is None:
            raise RuntimeError(f"Dataset not found: {job['dataset_id']}")

        data_path = dataset.get("path_h5ad") or dataset.get("path_raw")
        if not data_path:
            raise RuntimeError(f"No usable data path in dataset {dataset['id']}")

        params = {**job.get("params", {}), "data_path": data_path}
        control_dataset_id = params.get("control_dataset_id")
        if control_dataset_id:
            control_dataset = self.store.get_dataset(control_dataset_id)
            if control_dataset is None:
                raise RuntimeError(
                    f"Control dataset not found: {params['control_dataset_id']}"
                )
            params["control_data_path"] = control_dataset.get("path_h5ad") or control_dataset.get(
                "path_raw"
            )

        job_dir = settings.artifact_dir / "jobs" / job_id
        planned_log_path = job_dir / "train.log"
        self.store.update_job(
            job_id,
            {
                "artifact_dir": str(job_dir),
                "log_path": str(planned_log_path),
            },
        )

        def on_train_start(pid: int) -> None:
            self.store.update_job(job_id, {"train_pid": pid})

        train_out = self.runner.run_train(
            job_dir=job_dir,
            params=params,
            on_start=on_train_start,
        )
        model_record = self.store.create_model(
            {
                "job_id": job_id,
                "dataset_id": dataset["id"],
                "path_ckpt": train_out["model_path"],
                "log_path": train_out["log_path"],
                "gene_size": params["gene_size"],
                "output_dim": params["output_dim"],
                "use_drug_structure": bool(params.get("use_drug_structure", False)),
                "metrics": train_out.get("metrics", {}),
            }
        )
        self.store.update_job(
            job_id,
            {
                "model_id": model_record["id"],
                "artifact_dir": str(job_dir),
                "log_path": train_out["log_path"],
            },
        )

    def _execute_predict(self, job_id: str, job: dict[str, Any]) -> None:
        dataset = self.store.get_dataset(job["dataset_id"])
        if dataset is None:
            raise RuntimeError(f"Dataset not found: {job['dataset_id']}")

        data_path = dataset.get("path_h5ad") or dataset.get("path_raw")
        if not data_path:
            raise RuntimeError(f"No usable data path in dataset {dataset['id']}")

        params = {**job.get("params", {}), "data_path": data_path}
        model_id = params.get("model_id")
        if model_id and not params.get("model_path"):
            model = self.store.get_model(model_id)
            if model is None:
                raise RuntimeError(f"Model not found: {model_id}")
            params["model_path"] = model["path_ckpt"]

        if not params.get("model_path"):
            raise RuntimeError("model_path (or model_id) is required for predict job")

        job_dir = settings.artifact_dir / "jobs" / job_id
        pred_out = self.runner.run_predict(job_dir=job_dir, params=params)
        result_record = self.store.create_result(
            {
                "job_id": job_id,
                "dataset_id": dataset["id"],
                "prediction_path": pred_out["prediction_path"],
                "log_path": pred_out["log_path"],
                "summary": pred_out.get("summary", {}),
                "artifact_dir": str(job_dir),
            }
        )
        self.store.update_job(
            job_id,
            {
                "result_id": result_record["id"],
                "artifact_dir": str(job_dir),
                "log_path": pred_out["log_path"],
            },
        )
