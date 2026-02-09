from __future__ import annotations

from pathlib import Path

from fastapi.testclient import TestClient

from backend.app.api import seurat as seurat_api
from backend.app.main import app
from backend.app.runtime import store


def test_seurat_inspect_requires_existing_dataset() -> None:
    client = TestClient(app)
    response = client.post("/api/seurat/inspect", json={"dataset_id": "missing"})
    assert response.status_code == 404
    assert response.json()["detail"] == "Dataset not found"


def test_seurat_inspect_requires_h5ad_path() -> None:
    dataset = store.create_dataset(
        {
            "name": "no-h5ad",
            "status": "uploaded",
            "input_type": "rds",
            "path_raw": "fake.rds",
            "path_h5ad": None,
            "smiles_path": None,
            "validation": None,
        }
    )
    client = TestClient(app)
    response = client.post("/api/seurat/inspect", json={"dataset_id": dataset["id"]})
    assert response.status_code == 400
    assert "no prepared h5ad path" in response.json()["detail"]


def test_prepare_training_requires_existing_dataset() -> None:
    client = TestClient(app)
    response = client.post(
        "/api/seurat/prepare-training",
        json={
            "dataset_id": "missing",
            "group_column": "sample",
            "cluster_column": "celltype",
            "selected_clusters": ["A"],
            "seed": 42,
        },
    )
    assert response.status_code == 404
    assert response.json()["detail"] == "Dataset not found"


def test_prepare_training_requires_h5ad_path() -> None:
    dataset = store.create_dataset(
        {
            "name": "no-h5ad-prepare",
            "status": "uploaded",
            "input_type": "rds",
            "path_raw": "fake.rds",
            "path_h5ad": None,
            "smiles_path": None,
            "validation": None,
        }
    )
    client = TestClient(app)
    response = client.post(
        "/api/seurat/prepare-training",
        json={
            "dataset_id": dataset["id"],
            "group_column": "sample",
            "cluster_column": "celltype",
            "selected_clusters": ["A"],
            "seed": 42,
        },
    )
    assert response.status_code == 400
    assert "no prepared h5ad path" in response.json()["detail"]


def test_prepare_training_success_with_stubbed_preprocessor(
    monkeypatch,
    tmp_path: Path,
) -> None:
    source_h5ad = tmp_path / "source.h5ad"
    source_h5ad.write_text("placeholder", encoding="utf-8")

    dataset = store.create_dataset(
        {
            "name": "source",
            "status": "validated",
            "input_type": "h5ad",
            "path_raw": str(source_h5ad),
            "path_h5ad": str(source_h5ad),
            "smiles_path": None,
            "validation": {
                "valid": True,
                "errors": [],
                "warnings": [],
                "summary": {"n_cells": 1000, "n_genes": 2000, "obs_columns": []},
            },
        }
    )

    def _fake_prepare_training_dataset(
        input_path: Path,
        output_path: Path,
        group_column: str,
        cluster_column: str,
        selected_clusters: list[str],
        seed: int,
        max_cells: int,
        max_genes: int,
    ) -> dict[str, object]:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("prepared", encoding="utf-8")
        return {
            "path_h5ad": str(output_path),
            "n_cells": min(max_cells, 500),
            "n_genes": min(max_genes, 500),
            "sampling_report": {
                "mode": "stratified_sampling",
                "seed": seed,
                "input_cells": 1000,
                "output_cells": 500,
            },
            "gene_report": {
                "method": "wilcoxon_deg",
                "input_genes": 2000,
                "output_genes": 500,
                "selected_genes": ["G1", "G2"],
            },
            "group_column": group_column,
            "cluster_column": cluster_column,
            "selected_clusters": selected_clusters,
        }

    monkeypatch.setattr(
        seurat_api,
        "prepare_training_dataset",
        _fake_prepare_training_dataset,
    )

    client = TestClient(app)
    response = client.post(
        "/api/seurat/prepare-training",
        json={
            "dataset_id": dataset["id"],
            "group_column": "sample",
            "cluster_column": "celltype",
            "selected_clusters": ["T", "B"],
            "seed": 123,
        },
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["n_cells"] == 500
    assert payload["n_genes"] == 500
    assert payload["prepared_dataset_id"]
    assert payload["job_id"]

    job_response = client.get(f"/api/seurat/prepare-training/{payload['job_id']}")
    assert job_response.status_code == 200
    job = job_response.json()["job"]
    assert job["status"] == "success"
    assert job["prepared_dataset_id"] == payload["prepared_dataset_id"]
