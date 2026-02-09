from __future__ import annotations

from pathlib import Path
from typing import Any


def _load_adata(data_path: Path) -> Any:
    try:
        import scanpy as sc
    except ModuleNotFoundError as exc:  # noqa: PERF203
        raise RuntimeError(
            "scanpy is not installed in current environment. "
            "Install dependencies before Seurat inspection."
        ) from exc
    return sc.read_h5ad(str(data_path))


def _pick_umap_key(obsm: Any) -> str | None:
    keys = [str(key) for key in obsm.keys()]
    for candidate in ("X_umap", "X_UMAP", "umap", "UMAP"):
        if candidate in keys:
            return candidate
    for key in keys:
        if "umap" in key.lower():
            return key
    return None


def _embedding_shape(embedding: Any) -> tuple[int, int]:
    shape = getattr(embedding, "shape", None)
    if shape is not None and len(shape) >= 2:
        return int(shape[0]), int(shape[1])

    n_points = len(embedding)
    n_dims = len(embedding[0]) if n_points > 0 else 0
    return n_points, n_dims


def _embedding_xy(embedding: Any, row_index: int) -> tuple[float, float]:
    try:
        x_value = embedding[row_index, 0]
        y_value = embedding[row_index, 1]
    except Exception:  # noqa: BLE001
        row = embedding[row_index]
        x_value = row[0]
        y_value = row[1]

    x = float(x_value.item() if hasattr(x_value, "item") else x_value)
    y = float(y_value.item() if hasattr(y_value, "item") else y_value)
    return x, y


def _cell_id(obs: Any, row_index: int) -> str:
    index = getattr(obs, "index", None)
    if index is None:
        return str(row_index)
    try:
        return str(index[row_index])
    except Exception:  # noqa: BLE001
        return str(row_index)


def inspect_h5ad(data_path: Path, umap_preview_limit: int = 1500) -> dict[str, object]:
    if not data_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {data_path}")

    adata = _load_adata(data_path)
    obs = getattr(adata, "obs", None)

    metadata_columns: list[str] = []
    if obs is not None and hasattr(obs, "columns"):
        metadata_columns = sorted(str(name) for name in obs.columns.tolist())

    n_cells = int(getattr(adata, "n_obs", 0))
    n_genes = int(getattr(adata, "n_vars", 0))
    warnings: list[str] = []

    umap_payload: dict[str, object] | None = None
    obsm = getattr(adata, "obsm", None)
    if obsm is None or not hasattr(obsm, "keys"):
        warnings.append("UMAP embedding not found: AnnData.obsm is missing.")
    else:
        umap_key = _pick_umap_key(obsm)
        if umap_key is None:
            warnings.append("UMAP embedding not found in AnnData.obsm.")
        else:
            embedding = obsm[umap_key]
            n_points, n_dims = _embedding_shape(embedding)
            if n_dims < 2:
                warnings.append(
                    f"UMAP embedding '{umap_key}' has less than 2 dimensions."
                )
            else:
                preview_size = max(1, min(int(umap_preview_limit), n_points))
                preview: list[dict[str, object]] = []
                for row_index in range(preview_size):
                    x, y = _embedding_xy(embedding, row_index)
                    preview.append(
                        {
                            "cell_id": _cell_id(obs, row_index),
                            "x": x,
                            "y": y,
                        }
                    )

                umap_payload = {
                    "key": umap_key,
                    "n_points": n_points,
                    "preview_count": len(preview),
                    "truncated": n_points > len(preview),
                    "preview": preview,
                }

    return {
        "n_cells": n_cells,
        "n_genes": n_genes,
        "metadata_columns": metadata_columns,
        "has_umap": umap_payload is not None,
        "umap": umap_payload,
        "warnings": warnings,
    }
