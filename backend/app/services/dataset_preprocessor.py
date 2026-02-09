from __future__ import annotations

import math
import random
from pathlib import Path
from typing import Any

import numpy as np


def _load_adata(data_path: Path) -> Any:
    try:
        import scanpy as sc
    except ModuleNotFoundError as exc:  # noqa: PERF203
        raise RuntimeError(
            "scanpy is not installed in current environment. "
            "Install dependencies before prepare-training."
        ) from exc
    return sc.read_h5ad(str(data_path))


def _normalize_obs_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float) and math.isnan(value):
        return "NA"
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "<na>", "nat"}:
        return "NA"
    return text


def _obs_column_as_strings(obs: Any, column: str) -> list[str]:
    values = obs[column].tolist()
    return [_normalize_obs_value(item) for item in values]


def _count_values(values: list[str], indices: list[int]) -> dict[str, int]:
    counter: dict[str, int] = {}
    for idx in indices:
        key = values[idx]
        counter[key] = counter.get(key, 0) + 1
    return dict(sorted(counter.items(), key=lambda item: item[0]))


def _allocate_proportional(counts: dict[str, int], target: int) -> dict[str, int]:
    if target <= 0:
        return {key: 0 for key in counts}

    total = sum(counts.values())
    if total <= target:
        return dict(counts)

    allocated: dict[str, int] = {}
    remainders: list[tuple[float, int, str]] = []
    assigned = 0

    for key in sorted(counts.keys()):
        exact = counts[key] * target / total
        base = int(math.floor(exact))
        allocated[key] = min(base, counts[key])
        assigned += allocated[key]
        remainders.append((exact - base, counts[key] - allocated[key], key))

    remaining = target - assigned
    remainders.sort(key=lambda item: (-item[0], -item[1], item[2]))
    while remaining > 0:
        progressed = False
        for _, capacity, key in remainders:
            if remaining <= 0:
                break
            if capacity <= 0:
                continue
            if allocated[key] >= counts[key]:
                continue
            allocated[key] += 1
            remaining -= 1
            progressed = True
        if not progressed:
            break

    return allocated


def _sample_within_indices(
    rng: random.Random,
    indices: list[int],
    target: int,
) -> list[int]:
    if target >= len(indices):
        return list(indices)
    return rng.sample(indices, target)


def stratified_sample_cells(
    group_values: list[str],
    cluster_values: list[str],
    max_cells: int,
    seed: int,
) -> tuple[list[int], dict[str, object]]:
    all_indices = list(range(len(group_values)))
    if len(all_indices) <= max_cells:
        report = {
            "mode": "full",
            "seed": seed,
            "max_cells": max_cells,
            "input_cells": len(all_indices),
            "output_cells": len(all_indices),
            "group_counts_before": _count_values(group_values, all_indices),
            "group_counts_after": _count_values(group_values, all_indices),
            "cluster_counts_before": _count_values(cluster_values, all_indices),
            "cluster_counts_after": _count_values(cluster_values, all_indices),
        }
        return all_indices, report

    rng = random.Random(seed)
    target = max_cells
    group_to_indices: dict[str, list[int]] = {}
    for idx, group_name in enumerate(group_values):
        group_to_indices.setdefault(group_name, []).append(idx)

    group_counts = {group: len(indices) for group, indices in group_to_indices.items()}
    group_quota = _allocate_proportional(group_counts, target)

    selected: list[int] = []
    for group_name in sorted(group_to_indices.keys()):
        indices = group_to_indices[group_name]
        quota = group_quota.get(group_name, 0)
        if quota <= 0:
            continue
        if len(indices) <= quota:
            selected.extend(indices)
            continue

        cluster_to_indices: dict[str, list[int]] = {}
        for idx in indices:
            cluster_to_indices.setdefault(cluster_values[idx], []).append(idx)

        if len(cluster_to_indices) <= 1:
            selected.extend(_sample_within_indices(rng, indices, quota))
            continue

        cluster_counts = {
            cluster_name: len(cluster_indices)
            for cluster_name, cluster_indices in cluster_to_indices.items()
        }
        cluster_quota = _allocate_proportional(cluster_counts, quota)
        group_selected: list[int] = []
        for cluster_name in sorted(cluster_to_indices.keys()):
            cluster_indices = cluster_to_indices[cluster_name]
            cluster_target = cluster_quota.get(cluster_name, 0)
            if cluster_target <= 0:
                continue
            group_selected.extend(
                _sample_within_indices(rng, cluster_indices, cluster_target)
            )

        if len(group_selected) < quota:
            group_selected_set = set(group_selected)
            remaining = [idx for idx in indices if idx not in group_selected_set]
            group_selected.extend(
                _sample_within_indices(rng, remaining, quota - len(group_selected))
            )
        elif len(group_selected) > quota:
            group_selected = _sample_within_indices(rng, group_selected, quota)

        selected.extend(group_selected)

    if len(selected) < target:
        selected_set = set(selected)
        remaining = [idx for idx in all_indices if idx not in selected_set]
        selected.extend(_sample_within_indices(rng, remaining, target - len(selected)))
    elif len(selected) > target:
        selected = _sample_within_indices(rng, selected, target)

    selected = sorted(selected)
    report = {
        "mode": "stratified_sampling",
        "seed": seed,
        "max_cells": max_cells,
        "input_cells": len(all_indices),
        "output_cells": len(selected),
        "group_counts_before": _count_values(group_values, all_indices),
        "group_counts_after": _count_values(group_values, selected),
        "cluster_counts_before": _count_values(cluster_values, all_indices),
        "cluster_counts_after": _count_values(cluster_values, selected),
    }
    return selected, report


def _rank_genes_by_variance(adata: Any) -> list[str]:
    matrix = adata.X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix_np = np.asarray(matrix)
    if matrix_np.ndim == 1:
        matrix_np = matrix_np.reshape(1, -1)
    variances = np.var(matrix_np, axis=0)
    top_idx = np.argsort(-variances)
    return [str(adata.var_names[idx]) for idx in top_idx.tolist()]


def _fallback_hvg_genes(adata: Any, max_genes: int) -> tuple[list[str], str]:
    fallback_method = "variance"
    selected: list[str] = []

    try:
        import scanpy as sc

        work = adata.copy()
        sc.pp.highly_variable_genes(work, n_top_genes=max_genes, flavor="seurat")
        mask = work.var["highly_variable"] == True  # noqa: E712
        hvg = work.var[mask]
        if "highly_variable_rank" in hvg.columns:
            hvg = hvg.sort_values("highly_variable_rank")
        selected = [str(name) for name in hvg.index.tolist()[:max_genes]]
        if selected:
            fallback_method = "hvg"
            return selected, fallback_method
    except Exception:  # noqa: BLE001
        pass

    ranked = _rank_genes_by_variance(adata)
    selected = ranked[:max_genes]
    return selected, fallback_method


def _flatten_deg_candidates(
    rank_result: dict[str, Any],
) -> list[tuple[str, float, float]]:
    names = rank_result.get("names")
    pvals_adj = rank_result.get("pvals_adj")
    scores = rank_result.get("scores")
    if names is None:
        return []

    candidates: list[tuple[str, float, float]] = []
    dtype_names = getattr(getattr(names, "dtype", None), "names", None)
    if dtype_names:
        for group_name in dtype_names:
            group_genes = names[group_name]
            group_pvals = pvals_adj[group_name] if pvals_adj is not None else None
            group_scores = scores[group_name] if scores is not None else None
            for rank_idx, gene in enumerate(group_genes):
                gene_name = str(gene)
                pvalue = (
                    float(group_pvals[rank_idx])
                    if group_pvals is not None
                    else float("inf")
                )
                score = (
                    float(group_scores[rank_idx]) if group_scores is not None else 0.0
                )
                candidates.append((gene_name, pvalue, score))
        return candidates

    names_arr = np.asarray(names)
    pvals_arr = np.asarray(pvals_adj) if pvals_adj is not None else None
    scores_arr = np.asarray(scores) if scores is not None else None
    if names_arr.ndim == 1:
        for idx, gene in enumerate(names_arr.tolist()):
            pvalue = float(pvals_arr[idx]) if pvals_arr is not None else float("inf")
            score = float(scores_arr[idx]) if scores_arr is not None else 0.0
            candidates.append((str(gene), pvalue, score))
    else:
        n_rank, n_groups = names_arr.shape[0], names_arr.shape[1]
        for rank_idx in range(n_rank):
            for group_idx in range(n_groups):
                gene = str(names_arr[rank_idx, group_idx])
                pvalue = (
                    float(pvals_arr[rank_idx, group_idx])
                    if pvals_arr is not None
                    else float("inf")
                )
                score = (
                    float(scores_arr[rank_idx, group_idx])
                    if scores_arr is not None
                    else 0.0
                )
                candidates.append((gene, pvalue, score))

    return candidates


def select_top_genes(
    adata: Any,
    max_genes: int,
    seed: int,
) -> tuple[list[str], dict[str, object]]:
    n_genes = int(adata.n_vars)
    if n_genes <= max_genes:
        selected = [str(name) for name in adata.var_names.tolist()]
        return selected, {
            "method": "all_genes",
            "seed": seed,
            "input_genes": n_genes,
            "output_genes": len(selected),
            "fallback_used": False,
            "selected_genes": selected,
        }

    group_values = _obs_column_as_strings(adata.obs, "Group")
    group_counts = _count_values(group_values, list(range(len(group_values))))
    valid_deg_groups = len(group_counts) >= 2 and min(group_counts.values()) >= 2

    selected: list[str] = []
    fallback_used = False
    fallback_reason: str | None = None
    method = "wilcoxon_deg"

    if valid_deg_groups:
        try:
            import scanpy as sc

            rank_n_genes = min(n_genes, max_genes * 5)
            sc.tl.rank_genes_groups(
                adata,
                groupby="Group",
                method="wilcoxon",
                n_genes=rank_n_genes,
            )
            rank_result = adata.uns.get("rank_genes_groups", {})
            candidates = _flatten_deg_candidates(rank_result)
            candidates.sort(key=lambda item: (item[1], -abs(item[2]), item[0]))
            seen: set[str] = set()
            valid_names = {str(name) for name in adata.var_names.tolist()}
            for gene_name, _, _ in candidates:
                if gene_name not in valid_names or gene_name in seen:
                    continue
                seen.add(gene_name)
                selected.append(gene_name)
                if len(selected) >= max_genes:
                    break
            if len(selected) < max_genes:
                fallback_used = True
                fallback_reason = "deg_not_enough_genes"
        except Exception as exc:  # noqa: BLE001
            fallback_used = True
            fallback_reason = f"deg_failed: {exc}"
    else:
        fallback_used = True
        fallback_reason = "insufficient_groups_for_deg"

    if fallback_used:
        fallback_genes, fallback_method = _fallback_hvg_genes(adata, max_genes)
        method = fallback_method if not selected else f"wilcoxon_deg+{fallback_method}"
        selected_set = set(selected)
        for gene_name in fallback_genes:
            if gene_name in selected_set:
                continue
            selected.append(gene_name)
            selected_set.add(gene_name)
            if len(selected) >= max_genes:
                break

    selected = selected[:max_genes]
    return selected, {
        "method": method,
        "seed": seed,
        "input_genes": n_genes,
        "output_genes": len(selected),
        "group_counts": group_counts,
        "fallback_used": fallback_used,
        "fallback_reason": fallback_reason,
        "selected_genes": selected,
    }


def prepare_training_dataset(
    input_path: Path,
    output_path: Path,
    group_column: str,
    cluster_column: str,
    selected_clusters: list[str],
    seed: int = 42,
    max_cells: int = 500,
    max_genes: int = 500,
) -> dict[str, object]:
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")
    if not selected_clusters:
        raise ValueError("selected_clusters must not be empty")

    adata = _load_adata(input_path)
    obs_columns = set(adata.obs.columns.tolist())
    if group_column not in obs_columns:
        raise ValueError(f"group_column not found in obs: {group_column}")
    if cluster_column not in obs_columns:
        raise ValueError(f"cluster_column not found in obs: {cluster_column}")

    selected_cluster_values = {_normalize_obs_value(item) for item in selected_clusters}
    input_cluster_values = _obs_column_as_strings(adata.obs, cluster_column)
    filtered_indices = [
        idx
        for idx, cluster_name in enumerate(input_cluster_values)
        if cluster_name in selected_cluster_values
    ]
    if not filtered_indices:
        raise ValueError("No cells left after cluster filtering")

    adata_filtered = adata[filtered_indices].copy()
    group_values = _obs_column_as_strings(adata_filtered.obs, group_column)
    cluster_values = _obs_column_as_strings(adata_filtered.obs, cluster_column)
    adata_filtered.obs["Group"] = group_values
    adata_filtered.obs["Cluster"] = cluster_values

    sampled_indices, sampling_report = stratified_sample_cells(
        group_values=group_values,
        cluster_values=cluster_values,
        max_cells=max_cells,
        seed=seed,
    )
    adata_sampled = adata_filtered[sampled_indices].copy()

    selected_genes, gene_report = select_top_genes(
        adata=adata_sampled,
        max_genes=max_genes,
        seed=seed,
    )
    if not selected_genes:
        raise RuntimeError("No genes selected for training matrix")

    adata_final = adata_sampled[:, selected_genes].copy()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata_final.write_h5ad(str(output_path))

    return {
        "path_h5ad": str(output_path),
        "n_cells": int(adata_final.n_obs),
        "n_genes": int(adata_final.n_vars),
        "sampling_report": sampling_report,
        "gene_report": gene_report,
        "group_column": group_column,
        "cluster_column": cluster_column,
        "selected_clusters": sorted(selected_cluster_values),
    }
