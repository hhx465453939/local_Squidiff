from __future__ import annotations

from fastapi.testclient import TestClient

from backend.app.main import app
from backend.app.runtime import store


def _create_h5ad_dataset(
    *,
    name: str,
    path_suffix: str,
    prepared_from_dataset_id: str | None = None,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "name": name,
        "status": "prepared" if prepared_from_dataset_id else "validated",
        "input_type": "h5ad",
        "path_raw": f"tmp_cache/{path_suffix}.h5ad",
        "path_h5ad": f"tmp_cache/{path_suffix}.h5ad",
        "smiles_path": None,
        "validation": {
            "valid": True,
            "errors": [],
            "warnings": [],
            "summary": {"n_cells": 100, "n_genes": 200, "obs_columns": []},
        },
    }
    if prepared_from_dataset_id:
        payload["prepared_from_dataset_id"] = prepared_from_dataset_id
    return store.create_dataset(payload)


def test_train_job_auto_uses_latest_prepared_dataset() -> None:
    source = _create_h5ad_dataset(name="source-a", path_suffix="source-a")
    _create_h5ad_dataset(
        name="source-a-prepared-1",
        path_suffix="source-a-prepared-1",
        prepared_from_dataset_id=source["id"],
    )
    latest_prepared = _create_h5ad_dataset(
        name="source-a-prepared-2",
        path_suffix="source-a-prepared-2",
        prepared_from_dataset_id=source["id"],
    )

    client = TestClient(app)
    response = client.post(
        "/api/jobs/train",
        json={
            "dataset_id": source["id"],
            "gene_size": 100,
            "output_dim": 100,
            "batch_size": 32,
            "lr": 1e-4,
        },
    )
    assert response.status_code == 200
    job = response.json()["job"]
    assert job["dataset_id"] == latest_prepared["id"]
    assert job["source_dataset_id"] == source["id"]
    assert job["prepared_dataset_id"] == latest_prepared["id"]
    assert job["used_prepared_dataset"] is True
    assert job["params"]["requested_dataset_id"] == source["id"]
    assert job["params"]["train_dataset_id"] == latest_prepared["id"]
    assert job["params"]["auto_selected_prepared_dataset"] is True


def test_train_job_rejects_mismatched_prepared_dataset() -> None:
    source = _create_h5ad_dataset(name="source-b", path_suffix="source-b")
    other_source = _create_h5ad_dataset(name="source-c", path_suffix="source-c")
    other_prepared = _create_h5ad_dataset(
        name="source-c-prepared",
        path_suffix="source-c-prepared",
        prepared_from_dataset_id=other_source["id"],
    )

    client = TestClient(app)
    response = client.post(
        "/api/jobs/train",
        json={
            "dataset_id": source["id"],
            "prepared_dataset_id": other_prepared["id"],
            "gene_size": 100,
            "output_dim": 100,
            "batch_size": 32,
            "lr": 1e-4,
        },
    )
    assert response.status_code == 400
    assert "does not belong" in response.json()["detail"]
