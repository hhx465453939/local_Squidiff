from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

from ..runtime import store

router = APIRouter(prefix="/api/results", tags=["results"])


def _with_asset_urls(result: dict[str, object]) -> dict[str, object]:
    copy = dict(result)
    summary = copy.get("summary")
    if not isinstance(summary, dict):
        return copy
    assets = summary.get("assets")
    if not isinstance(assets, list):
        return copy

    mapped_assets: list[dict[str, object]] = []
    for item in assets:
        if not isinstance(item, dict):
            continue
        path = item.get("path")
        if not isinstance(path, str):
            continue
        name = Path(path).name
        mapped_assets.append(
            {
                **item,
                "name": name,
                "url": f"/api/results/{copy['id']}/assets/{name}",
            }
        )
    copy["summary"] = {**summary, "assets": mapped_assets}
    return copy


@router.get("")
async def list_results() -> dict[str, object]:
    return {"items": [_with_asset_urls(item) for item in store.list_results()]}


@router.get("/{result_id}")
async def get_result(result_id: str) -> dict[str, object]:
    result = store.get_result(result_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Result not found")
    return {"result": _with_asset_urls(result)}


@router.get("/job/{job_id}")
async def get_result_by_job(job_id: str) -> dict[str, object]:
    for result in store.list_results():
        if result.get("job_id") == job_id:
            return {"result": _with_asset_urls(result)}
    raise HTTPException(status_code=404, detail="Result not found for job")


@router.get("/models/list")
async def list_models() -> dict[str, object]:
    return {"items": store.list_models()}


@router.get("/models/{model_id}")
async def get_model(model_id: str) -> dict[str, object]:
    model = store.get_model(model_id)
    if model is None:
        raise HTTPException(status_code=404, detail="Model not found")
    return {"model": model}


def _validate_asset_path(result: dict[str, object], file_path: Path) -> Path:
    artifact_dir = result.get("artifact_dir")
    if not artifact_dir:
        raise HTTPException(status_code=404, detail="Result has no artifact directory")

    root = Path(str(artifact_dir)).resolve()
    target = file_path.resolve()
    if root not in target.parents and root != target:
        raise HTTPException(status_code=400, detail="Invalid asset path")
    if not target.exists():
        raise HTTPException(status_code=404, detail="Asset not found")
    return target


@router.get("/{result_id}/assets/{asset_name}")
async def get_result_asset(result_id: str, asset_name: str) -> FileResponse:
    result = store.get_result(result_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Result not found")

    summary = result.get("summary") or {}
    assets = summary.get("assets") if isinstance(summary, dict) else None
    if not isinstance(assets, list):
        raise HTTPException(status_code=404, detail="No assets found for result")

    for item in assets:
        if not isinstance(item, dict):
            continue
        path = item.get("path")
        if not isinstance(path, str):
            continue
        candidate = Path(path)
        if candidate.name != asset_name:
            continue
        safe_path = _validate_asset_path(result, candidate)
        return FileResponse(path=safe_path)

    raise HTTPException(status_code=404, detail="Asset not found")
