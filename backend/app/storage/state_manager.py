from __future__ import annotations

import json
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from uuid import uuid4


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


class JsonStateStore:
    """Simple JSON file store for MVP state persistence."""

    def __init__(self, state_dir: Path) -> None:
        self.state_dir = state_dir
        self._lock = threading.Lock()
        self._files = {
            "datasets": state_dir / "datasets.json",
            "jobs": state_dir / "jobs.json",
            "seurat_prepare_jobs": state_dir / "seurat_prepare_jobs.json",
            "models": state_dir / "models.json",
            "results": state_dir / "results.json",
        }
        self._ensure_files()

    def _ensure_files(self) -> None:
        self.state_dir.mkdir(parents=True, exist_ok=True)
        for file_path in self._files.values():
            if not file_path.exists():
                file_path.write_text("{}", encoding="utf-8")

    def _read_all(self, key: str) -> dict[str, Any]:
        file_path = self._files[key]
        raw = file_path.read_text(encoding="utf-8").strip()
        if not raw:
            return {}
        return json.loads(raw)

    def _write_all(self, key: str, data: dict[str, Any]) -> None:
        file_path = self._files[key]
        file_path.write_text(
            json.dumps(data, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

    def _create(self, key: str, payload: dict[str, Any]) -> dict[str, Any]:
        with self._lock:
            data = self._read_all(key)
            record_id = str(uuid4())
            record = {
                "id": record_id,
                "created_at": utc_now_iso(),
                "updated_at": utc_now_iso(),
                **payload,
            }
            data[record_id] = record
            self._write_all(key, data)
            return record

    def _get(self, key: str, record_id: str) -> dict[str, Any] | None:
        with self._lock:
            return self._read_all(key).get(record_id)

    def _list(self, key: str) -> list[dict[str, Any]]:
        with self._lock:
            return list(self._read_all(key).values())

    def _update(
        self,
        key: str,
        record_id: str,
        patch: dict[str, Any],
    ) -> dict[str, Any] | None:
        with self._lock:
            data = self._read_all(key)
            record = data.get(record_id)
            if record is None:
                return None
            record.update(patch)
            record["updated_at"] = utc_now_iso()
            data[record_id] = record
            self._write_all(key, data)
            return record

    def create_dataset(self, payload: dict[str, Any]) -> dict[str, Any]:
        return self._create("datasets", payload)

    def get_dataset(self, dataset_id: str) -> dict[str, Any] | None:
        return self._get("datasets", dataset_id)

    def list_datasets(self) -> list[dict[str, Any]]:
        return self._list("datasets")

    def update_dataset(
        self, dataset_id: str, patch: dict[str, Any]
    ) -> dict[str, Any] | None:
        return self._update("datasets", dataset_id, patch)

    def create_job(self, payload: dict[str, Any]) -> dict[str, Any]:
        return self._create("jobs", payload)

    def get_job(self, job_id: str) -> dict[str, Any] | None:
        return self._get("jobs", job_id)

    def list_jobs(self) -> list[dict[str, Any]]:
        return self._list("jobs")

    def update_job(self, job_id: str, patch: dict[str, Any]) -> dict[str, Any] | None:
        return self._update("jobs", job_id, patch)

    def create_seurat_prepare_job(self, payload: dict[str, Any]) -> dict[str, Any]:
        return self._create("seurat_prepare_jobs", payload)

    def get_seurat_prepare_job(self, job_id: str) -> dict[str, Any] | None:
        return self._get("seurat_prepare_jobs", job_id)

    def list_seurat_prepare_jobs(self) -> list[dict[str, Any]]:
        return self._list("seurat_prepare_jobs")

    def update_seurat_prepare_job(
        self,
        job_id: str,
        patch: dict[str, Any],
    ) -> dict[str, Any] | None:
        return self._update("seurat_prepare_jobs", job_id, patch)

    def create_model(self, payload: dict[str, Any]) -> dict[str, Any]:
        return self._create("models", payload)

    def get_model(self, model_id: str) -> dict[str, Any] | None:
        return self._get("models", model_id)

    def list_models(self) -> list[dict[str, Any]]:
        return self._list("models")

    def create_result(self, payload: dict[str, Any]) -> dict[str, Any]:
        return self._create("results", payload)

    def get_result(self, result_id: str) -> dict[str, Any] | None:
        return self._get("results", result_id)

    def list_results(self) -> list[dict[str, Any]]:
        return self._list("results")
