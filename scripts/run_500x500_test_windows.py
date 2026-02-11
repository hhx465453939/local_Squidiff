#!/usr/bin/env python3
"""
模拟前端 API 调用，用 data/ 下 TC.rds（筋膜）和 coTC.rds（结肠）测试 500×500 流程。
需先启动后端: uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000
Windows 下使用 conda r-4.3 进行 R 转换（cmd_conda）。
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

# 项目根目录
REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"

# Windows conda r-4.3（与 terminals 中 conda info --envs 一致）
DEFAULT_BASE_URL = "http://127.0.0.1:8000"
DEFAULT_R_CONDA_ENV = "r-4.3"
# Miniconda 在 F: 盘时的 conda.bat
DEFAULT_R_CONDA_BAT = os.getenv("LABFLOW_R_CONDA_BAT", "F:\\software\\Miniconda3\\condabin\\conda.bat")


def request_json(
    base_url: str,
    method: str,
    path: str,
    payload: dict | None = None,
    timeout_sec: int = 300,
) -> dict:
    url = f"{base_url.rstrip('/')}{path}"
    headers = {"Content-Type": "application/json"}
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    req = Request(url=url, data=data, headers=headers, method=method.upper())
    try:
        with urlopen(req, timeout=timeout_sec) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"{method} {path} -> {e.code}: {body}") from e
    except URLError as e:
        raise RuntimeError(f"{method} {path} -> {e.reason}") from e


def _infer_group_and_cluster_columns(h5ad_path: Path) -> tuple[str, str, list[str]]:
    """从 h5ad 的 obs 中推断 group_column、cluster_column 和 selected_clusters。"""
    try:
        import scanpy as sc
    except ModuleNotFoundError:
        raise RuntimeError("需要 scanpy 才能自动推断列名，请安装: pip install scanpy") from None
    adata = sc.read_h5ad(str(h5ad_path))
    cols = list(adata.obs.columns)
    # 常见命名
    group_candidates = ["Group", "group", "orig.ident", "condition", "sample"]
    cluster_candidates = ["seurat_clusters", "Cluster", "cluster", "celltype", "louvain"]
    group_column = next((c for c in group_candidates if c in cols), None)
    cluster_column = next((c for c in cluster_candidates if c in cols), None)
    if not cluster_column:
        # 取第一个看起来像分组的列（离散、取值数在 2~200）
        for c in cols:
            n = adata.obs[c].astype(str).nunique()
            if 2 <= n <= 200:
                cluster_column = c
                break
    if not cluster_column:
        cluster_column = cols[0] if cols else "seurat_clusters"
    if not group_column:
        group_column = cluster_column  # 退化为同一列
    values = adata.obs[cluster_column].astype(str).unique().tolist()
    selected = sorted(values)[:20]  # 最多选 20 个 cluster，避免过多
    return group_column, cluster_column, selected


def main() -> None:
    base_url = os.getenv("LABFLOW_BASE_URL", DEFAULT_BASE_URL)
    r_conda_env = os.getenv("LABFLOW_R_CONDA_ENV", DEFAULT_R_CONDA_ENV)
    r_conda_bat = os.getenv("LABFLOW_R_CONDA_BAT", DEFAULT_R_CONDA_BAT)

    datasets = [
        ("TC.rds", "TC-筋膜"),
        ("coTC.rds", "coTC-结肠"),
    ]

    print("=== 500×500 流程测试（模拟前端 API，Windows + r-4.3）===")
    print(f"BASE_URL={base_url}")
    print(f"R_CONDA_ENV={r_conda_env}")
    print(f"R_CONDA_BAT={r_conda_bat}")
    print()

    for filename, name in datasets:
        rds_path = DATA_DIR / filename
        if not rds_path.exists():
            print(f"[SKIP] {filename} 不存在于 {DATA_DIR}")
            continue

        print(f"--- {name} ({filename}) ---")

        # 1) 注册本地文件（复制到 backend/uploads）
        reg = request_json(
            base_url, "POST", "/api/datasets/register-local",
            payload={
                "local_path": str(rds_path),
                "dataset_name": name,
                "copy_to_upload_dir": True,
            },
            timeout_sec=60,
        )
        dataset_id = reg["dataset"]["id"]
        print(f"  注册: dataset_id={dataset_id}")

        # 2) 校验 + R 转 h5ad（cmd_conda + r-4.3）
        try:
            val = request_json(
                base_url, "POST", f"/api/datasets/{dataset_id}/validate",
                payload={
                    "auto_convert": True,
                    "r_exec_mode": "cmd_conda",
                    "r_conda_env": r_conda_env,
                    "r_conda_bat": r_conda_bat,
                },
                timeout_sec=600,
            )
        except Exception as e:
            print(f"  校验/转换失败: {e}")
            continue
        print(f"  校验: status={val['dataset'].get('status')}")

        path_h5ad = val["dataset"].get("path_h5ad")
        if not path_h5ad or not Path(path_h5ad).exists():
            print(f"  path_h5ad 缺失或文件不存在: {path_h5ad}")
            continue

        # 3) Inspect（可选，仅打印）
        try:
            insp = request_json(
                base_url, "POST", "/api/seurat/inspect",
                payload={"dataset_id": dataset_id, "umap_preview_limit": 500},
                timeout_sec=120,
            )
            meta = insp.get("inspect", {}).get("metadata_columns", [])
            print(f"  inspect: n_cells={insp.get('inspect', {}).get('n_cells')}, columns={meta[:8]}...")
        except Exception as e:
            print(f"  inspect 失败（继续）: {e}")

        # 4) 推断 group / cluster 列并调用 prepare-training 500×500
        try:
            group_col, cluster_col, selected = _infer_group_and_cluster_columns(Path(path_h5ad))
        except Exception as e:
            print(f"  推断列名失败: {e}")
            continue
        print(f"  使用: group_column={group_col}, cluster_column={cluster_col}, selected_clusters={selected[:5]}...")

        try:
            prep = request_json(
                base_url, "POST", "/api/seurat/prepare-training",
                payload={
                    "dataset_id": dataset_id,
                    "group_column": group_col,
                    "cluster_column": cluster_col,
                    "selected_clusters": selected,
                    "seed": 42,
                },
                timeout_sec=300,
            )
        except Exception as e:
            print(f"  prepare-training 失败: {e}")
            continue

        n_cells = prep.get("n_cells")
        n_genes = prep.get("n_genes")
        prep_id = prep.get("prepared_dataset_id")
        print(f"  prepare-training: prepared_dataset_id={prep_id}, n_cells={n_cells}, n_genes={n_genes}")
        if n_cells is not None and n_genes is not None:
            if 0 < n_cells <= 500 and 0 < n_genes <= 500:
                print("  -> 500×500 逻辑校验通过")
            else:
                print(f"  -> 注意: 未在 500×500 范围内 (cells={n_cells}, genes={n_genes})")
        print()

    print("=== 测试结束 ===")


if __name__ == "__main__":
    main()
    sys.exit(0)
