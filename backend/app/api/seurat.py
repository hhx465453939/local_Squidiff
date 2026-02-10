from __future__ import annotations

from pathlib import Path
from uuid import uuid4

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from ..core.config import settings
from ..runtime import store
from ..services.dataset_preprocessor import prepare_training_dataset
from ..services.seurat_inspector import inspect_h5ad
from ..storage.state_manager import utc_now_iso

router = APIRouter(prefix="/api/seurat", tags=["seurat"])


class SeuratInspectPayload(BaseModel):
    dataset_id: str
    umap_preview_limit: int = Field(default=1500, ge=1, le=5000)


class SeuratPrepareTrainingPayload(BaseModel):
    dataset_id: str
    group_column: str
    cluster_column: str
    selected_clusters: list[str] = Field(min_length=1)
    seed: int = Field(default=42, ge=0)


def _require_dataset(dataset_id: str) -> dict[str, object]:
    dataset = store.get_dataset(dataset_id)
    if dataset is None:
        raise HTTPException(status_code=404, detail="Dataset not found")
    return dataset


def _require_h5ad_path(dataset: dict[str, object]) -> Path:
    h5ad_path_raw = dataset.get("path_h5ad")
    if not isinstance(h5ad_path_raw, str) or not h5ad_path_raw:
        raise HTTPException(
            status_code=400,
            detail=(
                "Dataset has no prepared h5ad path. "
                "Please run /api/datasets/{dataset_id}/validate with auto_convert first."
            ),
        )
    return Path(h5ad_path_raw)


@router.post("/inspect")
async def inspect_seurat(payload: SeuratInspectPayload) -> dict[str, object]:
    dataset = _require_dataset(payload.dataset_id)
    h5ad_path = _require_h5ad_path(dataset)

    try:
        report = inspect_h5ad(
            data_path=h5ad_path,
            umap_preview_limit=payload.umap_preview_limit,
        )
    except FileNotFoundError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except RuntimeError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    return {
        "dataset_id": payload.dataset_id,
        "inspect": report,
    }


@router.post("/prepare-training")
async def prepare_training(
    payload: SeuratPrepareTrainingPayload,
) -> dict[str, object]:
    dataset = _require_dataset(payload.dataset_id)
    h5ad_path = _require_h5ad_path(dataset)

    job = store.create_seurat_prepare_job(
        {
            "status": "running",
            "dataset_id": payload.dataset_id,
            "params": payload.model_dump(),
            "error_msg": None,
            "started_at": utc_now_iso(),
            "ended_at": None,
            "prepared_dataset_id": None,
            "result": None,
        }
    )

    output_dir = settings.upload_dir / "prepared"
    output_path = (
        output_dir / f"{payload.dataset_id}-{job['id'][:8]}-{uuid4().hex[:8]}.h5ad"
    )

    try:
        result = prepare_training_dataset(
            input_path=h5ad_path,
            output_path=output_path,
            group_column=payload.group_column,
            cluster_column=payload.cluster_column,
            selected_clusters=payload.selected_clusters,
            seed=payload.seed,
            max_cells=500,
            max_genes=500,
        )

        prepared_dataset = store.create_dataset(
            {
                "name": f"{dataset.get('name', payload.dataset_id)}-prepared",
                "status": "prepared",
                "input_type": "h5ad",
                "path_raw": result["path_h5ad"],
                "path_h5ad": result["path_h5ad"],
                "smiles_path": dataset.get("smiles_path"),
                "validation": {
                    "valid": True,
                    "errors": [],
                    "warnings": [],
                    "summary": {
                        "n_cells": result["n_cells"],
                        "n_genes": result["n_genes"],
                        "obs_columns": ["Group", "Cluster"],
                    },
                },
                "prepared_from_dataset_id": payload.dataset_id,
                "prepare_job_id": job["id"],
                "prepare_report": {
                    "sampling_report": result["sampling_report"],
                    "gene_report": result["gene_report"],
                },
            }
        )

        final_patch = {
            "status": "success",
            "prepared_dataset_id": prepared_dataset["id"],
            "result": {
                "prepared_dataset_id": prepared_dataset["id"],
                "n_cells": result["n_cells"],
                "n_genes": result["n_genes"],
                "sampling_report": result["sampling_report"],
                "gene_report": result["gene_report"],
            },
            "ended_at": utc_now_iso(),
        }
        job = store.update_seurat_prepare_job(job["id"], final_patch) or {
            **job,
            **final_patch,
        }
    except Exception as exc:  # noqa: BLE001
        final_patch = {
            "status": "failed",
            "error_msg": str(exc),
            "ended_at": utc_now_iso(),
        }
        store.update_seurat_prepare_job(job["id"], final_patch)
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    result_payload = job.get("result") if isinstance(job, dict) else None
    if not isinstance(result_payload, dict):
        raise HTTPException(
            status_code=500, detail="prepare-training job result missing"
        )

    return {
        "job_id": job["id"],
        **result_payload,
    }


@router.get("/prepare-training/{job_id}")
async def get_prepare_training_job(job_id: str) -> dict[str, object]:
    job = store.get_seurat_prepare_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Prepare-training job not found")
    return {"job": job}
