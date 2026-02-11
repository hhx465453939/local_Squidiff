from __future__ import annotations

import os
import signal
import subprocess
import traceback
from datetime import datetime, timezone
from queue import Queue
from threading import Event, Thread
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


class JobQueue:
    def __init__(self, store: JsonStateStore, runner: Any) -> None:
        self.store = store
        self.runner = runner
        self._queue: Queue[str] = Queue()
        self._stop_event = Event()
        self._worker: Thread | None = None

    def start(self) -> None:
        if self._worker and self._worker.is_alive():
            return
        self._stop_event.clear()
        self._worker = Thread(target=self._work_loop, daemon=True)
        self._worker.start()

    def stop(self) -> None:
        self._stop_event.set()
        self._queue.put("__STOP__")
        if self._worker:
            self._worker.join(timeout=3)

    def enqueue(self, job_id: str) -> None:
        self._queue.put(job_id)

    def cancel(self, job_id: str) -> dict[str, Any] | None:
        job = self.store.get_job(job_id)
        if job is None:
            return None

        status = str(job.get("status") or "")
        if status in {"success", "failed", "canceled"}:
            return job

        if status == "queued":
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

    def _work_loop(self) -> None:
        while not self._stop_event.is_set():
            job_id = self._queue.get()
            if job_id == "__STOP__":
                self._queue.task_done()
                continue

            job = self.store.get_job(job_id)
            if job is None:
                self._queue.task_done()
                continue

            if job.get("status") == "canceled":
                self._queue.task_done()
                continue

            try:
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
                    self.store.update_job(
                        job_id, {"status": "success", "ended_at": _now()}
                    )
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
                self._queue.task_done()

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
            params["control_data_path"] = control_dataset.get(
                "path_h5ad"
            ) or control_dataset.get("path_raw")

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
