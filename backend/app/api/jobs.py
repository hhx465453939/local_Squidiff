from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from ..runtime import job_queue, store

router = APIRouter(prefix="/api/jobs", tags=["jobs"])


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

    log_path = job.get("log_path")
    if not isinstance(log_path, str):
        raise HTTPException(status_code=404, detail="Log not found")

    path = Path(log_path)
    if not path.exists():
        raise HTTPException(status_code=404, detail="Log file not found")

    text = path.read_text(encoding="utf-8", errors="replace")
    if len(text) > 20000:
        text = text[-20000:]
    return {"log": text}


@router.post("/train")
async def submit_train_job(payload: TrainJobPayload) -> dict[str, object]:
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
async def submit_predict_job(payload: PredictJobPayload) -> dict[str, object]:
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
            "params": payload.model_dump(),
            "error_msg": None,
            "started_at": None,
            "ended_at": None,
        }
    )
    job_queue.enqueue(job["id"])
    return {"job": job}
