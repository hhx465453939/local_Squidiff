from __future__ import annotations

from datetime import datetime, timezone
import os
import shutil
import subprocess
from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel, Field

from ..auth import require_auth
from ..core.config import settings
from ..runtime import job_queue, store

router = APIRouter(prefix="/api/jobs", tags=["jobs"])


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class TrainJobPayload(BaseModel):
    dataset_id: str
    prepared_dataset_id: str | None = None
    gene_size: int = Field(gt=0)
    output_dim: int = Field(gt=0)
    use_drug_structure: bool = False
    control_dataset_id: str | None = None
    batch_size: int = 64
    lr: float = 1e-4


class PredictJobPayload(BaseModel):
    dataset_id: str
    model_id: str | None = None
    model_path: str | None = None
    gene_size: int = Field(gt=0)
    output_dim: int = Field(gt=0)
    use_drug_structure: bool = False


def _latest_prepared_dataset(source_dataset_id: str) -> dict[str, object] | None:
    candidates = [
        item
        for item in store.list_datasets()
        if item.get("prepared_from_dataset_id") == source_dataset_id
        and isinstance(item.get("path_h5ad"), str)
        and item.get("path_h5ad")
    ]
    if not candidates:
        return None
    candidates.sort(
        key=lambda item: str(item.get("updated_at") or item.get("created_at") or ""),
        reverse=True,
    )
    return candidates[0]


@router.get("")
async def list_jobs() -> dict[str, object]:
    return {"items": store.list_jobs()}


def _remove_tree_if_allowed(path: Path) -> None:
    try:
        root = settings.artifact_dir.resolve()
        target = path.resolve()
    except OSError:
        return
    if root not in target.parents and root != target:
        return
    if target.exists() and target.is_dir():
        shutil.rmtree(target, ignore_errors=True)


def _is_pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        if os.name == "nt":
            result = subprocess.run(
                ["tasklist", "/FI", f"PID eq {pid}", "/FO", "CSV", "/NH"],
                capture_output=True,
                text=True,
                timeout=6,
                check=False,
            )
            output = (result.stdout or "").strip()
            if not output or "No tasks are running" in output:
                return False
            return str(pid) in output
        os.kill(pid, 0)
        return True
    except (OSError, subprocess.SubprocessError):
        return False


def _is_job_actually_running(job: dict[str, object]) -> bool:
    if job.get("status") != "running":
        return False
    pid = job.get("train_pid")
    return isinstance(pid, int) and _is_pid_alive(pid)


def _cascade_delete_job(job_id: str, *, purge_artifacts: bool) -> dict[str, int]:
    related_model_ids = [
        str(model.get("id"))
        for model in store.list_models()
        if model.get("job_id") == job_id and model.get("id")
    ]
    related_result_ids = [
        str(result.get("id"))
        for result in store.list_results()
        if result.get("job_id") == job_id and result.get("id")
    ]

    for model_id in related_model_ids:
        store.delete_model(model_id)
    for result_id in related_result_ids:
        store.delete_result(result_id)
    store.delete_job(job_id)

    if purge_artifacts:
        _remove_tree_if_allowed(settings.artifact_dir / "jobs" / job_id)

    return {
        "removed_models": len(related_model_ids),
        "removed_results": len(related_result_ids),
    }


@router.get("/{job_id}")
async def get_job(job_id: str) -> dict[str, object]:
    job = store.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return {"job": job}


@router.get("/{job_id}/log")
async def get_job_log(job_id: str) -> dict[str, str]:
    job = store.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    def _read_tail(path: Path) -> str:
        if not path.exists() or not path.is_file():
            return ""
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return ""
        if len(text) > 20000:
            text = text[-20000:]
        return text

    candidates: list[Path] = []
    log_path = job.get("log_path")
    if isinstance(log_path, str) and log_path.strip():
        primary = Path(log_path)
        candidates.append(primary)
        candidates.append(primary.parent / "logger" / "log.txt")
    else:
        primary = None

    candidates.append(settings.artifact_dir / "jobs" / job_id / "logger" / "log.txt")
    candidates.append(settings.artifact_dir / "jobs" / job_id / "train.log")

    seen: set[str] = set()
    existing: list[Path] = []
    for candidate in candidates:
        key = str(candidate.resolve()) if candidate.exists() else str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.exists():
            existing.append(candidate)

    if not existing:
        return {"log": ""}

    primary_text = _read_tail(existing[0])
    if primary_text.strip():
        return {"log": primary_text}

    for path in existing[1:]:
        fallback_text = _read_tail(path)
        if fallback_text.strip():
            return {"log": fallback_text}

    return {"log": primary_text}


@router.post("/{job_id}/cancel")
async def cancel_job(job_id: str) -> dict[str, object]:
    job = store.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.get("status") in {"success", "failed", "canceled"}:
        raise HTTPException(
            status_code=400,
            detail=f"Job is already finished: {job.get('status')}",
        )

    if job.get("status") == "running" and not _is_job_actually_running(job):
        stale = store.update_job(
            job_id,
            {
                "status": "canceled",
                "ended_at": _now(),
                "error_msg": "Canceled stale job (process not found).",
                "cancel_requested": True,
                "cancel_signal_sent": False,
            },
        )
        if stale is None:
            raise HTTPException(status_code=404, detail="Job not found")
        return {"job": stale}

    updated = job_queue.cancel(job_id)
    if updated is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return {"job": updated}


@router.delete("/{job_id}")
async def delete_job(
    job_id: str,
    purge_artifacts: bool = Query(
        True,
        description="If true, also remove related files under backend/artifacts/jobs/{job_id}",
    ),
) -> dict[str, object]:
    job = store.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    status = str(job.get("status") or "")
    if status in {"queued", "running"}:
        if _is_job_actually_running(job):
            raise HTTPException(
                status_code=400,
                detail=(
                    "Job is still active. Stop/cancel it first, or use /api/jobs/flush "
                    "to clear stale tasks."
                ),
            )
        store.update_job(
            job_id,
            {
                "status": "canceled",
                "ended_at": _now(),
                "error_msg": "Deleted as stale active job.",
                "cancel_requested": True,
            },
        )

    removed = _cascade_delete_job(job_id, purge_artifacts=purge_artifacts)

    return {
        "deleted": True,
        "removed_models": removed["removed_models"],
        "removed_results": removed["removed_results"],
    }


@router.post("/flush")
async def flush_jobs(
    scope: str = Query(
        "active",
        pattern="^(active|all)$",
        description="active: queued/running only; all: delete all jobs",
    ),
    purge_artifacts: bool = Query(
        True,
        description="If true, remove backend/artifacts/jobs/{job_id} while clearing",
    ),
    force: bool = Query(
        True,
        description="If true, attempt to kill running process tree before clearing",
    ),
) -> dict[str, object]:
    jobs = store.list_jobs()
    if scope == "active":
        target_jobs = [
            item for item in jobs if item.get("status") in {"queued", "running"}
        ]
    else:
        target_jobs = jobs

    deleted_jobs = 0
    removed_models = 0
    removed_results = 0
    skipped_running: list[str] = []

    for job in target_jobs:
        job_id = str(job.get("id") or "")
        if not job_id:
            continue
        latest = store.get_job(job_id)
        if latest is None:
            continue

        status = str(latest.get("status") or "")
        if status in {"queued", "running"}:
            if force:
                job_queue.cancel(job_id)
                latest = store.get_job(job_id) or latest
            if _is_job_actually_running(latest):
                skipped_running.append(job_id)
                continue
            store.update_job(
                job_id,
                {
                    "status": "canceled",
                    "ended_at": _now(),
                    "error_msg": "Cleared by one-click flush.",
                    "cancel_requested": True,
                },
            )

        removed = _cascade_delete_job(job_id, purge_artifacts=purge_artifacts)
        deleted_jobs += 1
        removed_models += removed["removed_models"]
        removed_results += removed["removed_results"]

    return {
        "scope": scope,
        "deleted_jobs": deleted_jobs,
        "removed_models": removed_models,
        "removed_results": removed_results,
        "skipped_running_jobs": skipped_running,
    }


@router.post("/train")
async def submit_train_job(
    payload: TrainJobPayload,
    current_user: dict[str, object] | None = Depends(require_auth),
) -> dict[str, object]:
    dataset = store.get_dataset(payload.dataset_id)
    if dataset is None:
        raise HTTPException(status_code=404, detail="Dataset not found")

    requested_dataset_id = payload.dataset_id
    train_dataset = dataset
    source_dataset_id = requested_dataset_id
    prepared_dataset_id: str | None = None
    used_prepared_dataset = False

    requested_from_id = dataset.get("prepared_from_dataset_id")
    if isinstance(requested_from_id, str) and requested_from_id:
        source_dataset_id = requested_from_id
        prepared_dataset_id = dataset["id"]
        used_prepared_dataset = True
    elif payload.prepared_dataset_id:
        prepared_dataset = store.get_dataset(payload.prepared_dataset_id)
        if prepared_dataset is None:
            raise HTTPException(status_code=404, detail="Prepared dataset not found")
        if prepared_dataset.get("prepared_from_dataset_id") != requested_dataset_id:
            raise HTTPException(
                status_code=400,
                detail=(
                    "prepared_dataset_id does not belong to the provided dataset_id. "
                    "Please use a matching prepared dataset."
                ),
            )
        train_dataset = prepared_dataset
        prepared_dataset_id = prepared_dataset["id"]
        used_prepared_dataset = True
    else:
        auto_prepared_dataset = _latest_prepared_dataset(requested_dataset_id)
        if auto_prepared_dataset is not None:
            train_dataset = auto_prepared_dataset
            prepared_dataset_id = auto_prepared_dataset["id"]
            used_prepared_dataset = True

    params = payload.model_dump(exclude_none=True)
    params["requested_dataset_id"] = requested_dataset_id
    params["train_dataset_id"] = train_dataset["id"]
    params["auto_selected_prepared_dataset"] = (
        payload.prepared_dataset_id is None and used_prepared_dataset
    )

    job = store.create_job(
        {
            "type": "train",
            "status": "queued",
            "dataset_id": train_dataset["id"],
            "owner_user_id": int(current_user["id"]) if current_user else None,
            "source_dataset_id": source_dataset_id,
            "prepared_dataset_id": prepared_dataset_id,
            "used_prepared_dataset": used_prepared_dataset,
            "params": params,
            "error_msg": None,
            "started_at": None,
            "ended_at": None,
        }
    )
    job_queue.enqueue(job["id"])
    return {"job": job}


@router.post("/predict")
async def submit_predict_job(
    payload: PredictJobPayload,
    current_user: dict[str, object] | None = Depends(require_auth),
) -> dict[str, object]:
    dataset = store.get_dataset(payload.dataset_id)
    if dataset is None:
        raise HTTPException(status_code=404, detail="Dataset not found")

    if payload.model_id is None and payload.model_path is None:
        raise HTTPException(
            status_code=400, detail="model_id or model_path is required"
        )

    job = store.create_job(
        {
            "type": "predict",
            "status": "queued",
            "dataset_id": payload.dataset_id,
            "owner_user_id": int(current_user["id"]) if current_user else None,
            "params": payload.model_dump(),
            "error_msg": None,
            "started_at": None,
            "ended_at": None,
        }
    )
    job_queue.enqueue(job["id"])
    return {"job": job}
