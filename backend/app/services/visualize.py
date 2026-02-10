from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def _safe_matrix(x: np.ndarray) -> np.ndarray:
    if x.ndim == 1:
        return x.reshape(-1, 1)
    return x


def generate_visual_assets(
    *,
    prediction_matrix: np.ndarray,
    output_dir: Path,
) -> dict[str, object]:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "matplotlib is required to generate visualization assets. "
            "Please install matplotlib or disable visualization generation."
        ) from exc

    try:
        from sklearn.decomposition import PCA
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "scikit-learn is required to generate visualization assets. "
            "Please install scikit-learn or disable visualization generation."
        ) from exc

    output_dir.mkdir(parents=True, exist_ok=True)
    matrix = _safe_matrix(prediction_matrix)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(matrix)
    pca_path = output_dir / "pca_scatter.png"
    plt.figure(figsize=(8, 6))
    plt.scatter(coords[:, 0], coords[:, 1], s=6, alpha=0.6)
    plt.title("Prediction PCA")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig(pca_path, dpi=160)
    plt.close()

    top_k = min(30, matrix.shape[1])
    top_idx = np.argsort(np.var(matrix, axis=0))[-top_k:]
    heat = matrix[:, top_idx]
    heatmap_path = output_dir / "heatmap_top_var_genes.png"
    plt.figure(figsize=(10, 5))
    plt.imshow(heat, aspect="auto", interpolation="nearest")
    plt.title("Prediction Heatmap (Top Variance Genes)")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(heatmap_path, dpi=160)
    plt.close()

    summary = {
        "n_cells": int(matrix.shape[0]),
        "n_genes": int(matrix.shape[1]),
        "mean_expression": float(np.mean(matrix)),
        "std_expression": float(np.std(matrix)),
        "assets": [
            {"type": "pca", "path": str(pca_path)},
            {"type": "heatmap", "path": str(heatmap_path)},
        ],
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
