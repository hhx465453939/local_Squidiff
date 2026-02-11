#!/usr/bin/env python3
"""
全流程测试：RDS → h5ad 转换 → 500×500 预处理 → 训练 → 预测 → 可视化报告。
模拟前端 API 调用，从转换到出图全部走通。metadata 已写入 prepared h5ad 的 .obs（Group/Cluster），
训练使用该 h5ad 的表达矩阵（train_squidiff 不单独读 metadata，数据由 prepare-training 已定）。
需先启动后端。
- 若仅验证流程：后端设置 LABFLOW_DRY_RUN=true，则训练不执行（只写占位 model.pt），预测用随机矩阵，但可视化图会正常生成。
- 若需真实训练与真实模型：后端不要设置 LABFLOW_DRY_RUN（或设为 false），再跑本脚本；训练会调用 train_squidiff.py，模型落在 backend/artifacts/jobs/<train_job_id>/checkpoints/ 下。

- 训练轮询：首轮等待 TRAIN_POLL_TIMEOUT_SEC（默认 3600s）；若超时则判断是否延长：后端若返回本任务的 train_pid，
  则仅当该 PID 仍在 nvidia-smi 的 GPU 进程列表中时才延长；否则退化为「GPU 利用率高于阈值或任意 Python 进程在 GPU」。
  默认每次延长 30 分钟，总时长上限 4 小时。需本机有 nvidia-smi。环境变量：LABFLOW_TRAIN_GPU_BUSY_THRESHOLD、LABFLOW_TRAIN_EXTEND_SEC、LABFLOW_TRAIN_MAX_TOTAL_SEC。

故障排除：
- 端口被占用：若 8000 已被占用，可改在其它端口启动后端（如 8002），
  并设置 LABFLOW_BASE_URL=http://127.0.0.1:8002 再运行本脚本。
- validate 报 conda.bat 无法识别：说明当前监听的后端是旧进程，未包含 R 转换的临时 .bat 修复。
  请关闭旧后端后重新启动（或在新端口启动新后端），再运行本脚本。
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
OUTPUT_REPORT_DIR = REPO_ROOT / "scripts" / "output" / "full_flow_report"

DEFAULT_BASE_URL = "http://127.0.0.1:8000"
DEFAULT_R_CONDA_ENV = "r-4.3"
DEFAULT_R_CONDA_BAT = os.getenv(
    "LABFLOW_R_CONDA_BAT", "F:\\software\\Miniconda3\\condabin\\conda.bat"
)
POLL_INTERVAL_SEC = 5
TRAIN_POLL_TIMEOUT_SEC = 3600
PREDICT_POLL_TIMEOUT_SEC = 600
# 训练超时后若 GPU 仍忙则自动延长的参数
TRAIN_GPU_BUSY_THRESHOLD = int(os.getenv("LABFLOW_TRAIN_GPU_BUSY_THRESHOLD", "10"))
TRAIN_EXTEND_SEC = int(os.getenv("LABFLOW_TRAIN_EXTEND_SEC", "1800"))
TRAIN_MAX_TOTAL_SEC = int(os.getenv("LABFLOW_TRAIN_MAX_TOTAL_SEC", "14400"))


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
    try:
        import scanpy as sc
    except ModuleNotFoundError:
        raise RuntimeError(
            "需要 scanpy 才能自动推断列名，请安装: pip install scanpy"
        ) from None
    adata = sc.read_h5ad(str(h5ad_path))
    cols = list(adata.obs.columns)
    group_candidates = ["Group", "group", "orig.ident", "condition", "sample"]
    cluster_candidates = [
        "seurat_clusters",
        "Cluster",
        "cluster",
        "celltype",
        "louvain",
    ]
    group_column = next((c for c in group_candidates if c in cols), None)
    cluster_column = next((c for c in cluster_candidates if c in cols), None)
    if not cluster_column:
        for c in cols:
            n = adata.obs[c].astype(str).nunique()
            if 2 <= n <= 200:
                cluster_column = c
                break
    if not cluster_column:
        cluster_column = cols[0] if cols else "seurat_clusters"
    if not group_column:
        group_column = cluster_column
    values = adata.obs[cluster_column].astype(str).unique().tolist()
    selected = sorted(values)[:20]
    return group_column, cluster_column, selected


def get_gpu_utilization() -> float | None:
    """当前 GPU 利用率 0–100，无法获取时返回 None（无 nvidia-smi 或非 NVIDIA）。"""
    try:
        out = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=utilization.gpu",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if out.returncode != 0 or not out.stdout:
            return None
        line = out.stdout.strip().split("\n")[0].strip()
        match = re.search(r"\d+", line)
        return float(match.group(0)) if match else None
    except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
        return None


def get_gpu_process_names() -> list[str]:
    """当前占用 GPU 的进程名列表（用于判断是否有 Python/训练进程在跑）。"""
    names: list[str] = []
    try:
        out = subprocess.run(
            [
                "nvidia-smi",
                "--query-compute-apps=process_name",
                "--format=csv,noheader",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if out.returncode == 0 and out.stdout:
            for line in out.stdout.strip().split("\n"):
                name = line.strip().strip('"').strip()
                if name:
                    names.append(name)
            if names:
                return names
        # 部分驱动/WDDM 下 --query-compute-apps 可能为空，用 -q 输出解析
        qout = subprocess.run(
            ["nvidia-smi", "-q"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if qout.returncode != 0 or not qout.stdout:
            return names
        # 解析 "Process name" 行
        for line in qout.stdout.split("\n"):
            if "Process name" in line:
                parts = line.split(":", 1)
                if len(parts) == 2 and parts[1].strip():
                    names.append(parts[1].strip())
        return names
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return []


def get_gpu_pids() -> list[int]:
    """当前占用 GPU 的进程 PID 列表（用于精确判断本任务训练进程是否仍在跑）。"""
    pids: list[int] = []
    try:
        out = subprocess.run(
            [
                "nvidia-smi",
                "--query-compute-apps=pid",
                "--format=csv,noheader",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if out.returncode == 0 and out.stdout:
            for line in out.stdout.strip().split("\n"):
                s = line.strip().strip('"').strip()
                if s.isdigit():
                    pids.append(int(s))
            if pids:
                return pids
        qout = subprocess.run(
            ["nvidia-smi", "-q"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if qout.returncode != 0 or not qout.stdout:
            return pids
        for line in qout.stdout.split("\n"):
            if "Process ID" in line:
                parts = line.split(":", 1)
                if len(parts) == 2:
                    s = parts[1].strip()
                    if s.isdigit():
                        pids.append(int(s))
        return pids
    except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
        return []


def is_python_process_on_gpu() -> bool:
    """GPU 上是否有名称含 python 的进程（训练多为 python 进程）。"""
    for name in get_gpu_process_names():
        if "python" in name.lower():
            return True
    return False


def poll_job_until_done(
    base_url: str,
    job_id: str,
    poll_interval_sec: int,
    max_wait_sec: int,
    timeout_sec: int = 60,
    *,
    extend_on_gpu_busy: bool = False,
    gpu_busy_threshold: int = 10,
    extend_sec: int = 1800,
    max_total_sec: int = 14400,
) -> dict:
    deadline = time.time() + max_wait_sec
    started = time.time()
    status = "running"
    while True:
        if time.time() > deadline:
            if not extend_on_gpu_busy:
                payload = request_json(
                    base_url, "GET", f"/api/jobs/{job_id}", timeout_sec=timeout_sec
                )
                job = payload.get("job", {})
                status = str(job.get("status", ""))
                if status in ("success", "failed"):
                    return job
                raise RuntimeError(
                    f"Job {job_id} 未在 {max_wait_sec}s 内完成，当前 status={status}"
                )
            if time.time() - started >= max_total_sec:
                payload = request_json(
                    base_url, "GET", f"/api/jobs/{job_id}", timeout_sec=timeout_sec
                )
                job = payload.get("job", {})
                status = str(job.get("status", ""))
                if status in ("success", "failed"):
                    return job
                raise RuntimeError(
                    f"Job {job_id} 在 {max_total_sec}s 总时长内未完成，当前 status={status}"
                )
            # 先拉取 job，若有 train_pid 则仅当该 PID 仍在 GPU 上时才延长
            payload = request_json(
                base_url, "GET", f"/api/jobs/{job_id}", timeout_sec=timeout_sec
            )
            job = payload.get("job", {})
            status = str(job.get("status", ""))
            if status in ("success", "failed"):
                return job
            train_pid_raw = job.get("train_pid")
            train_pid_val: int | None = None
            if train_pid_raw is not None:
                try:
                    train_pid_val = int(train_pid_raw)
                except (TypeError, ValueError):
                    pass
            gpu_pids = get_gpu_pids()
            util: float | None = None
            python_on_gpu = False
            if train_pid_val is not None:
                extend = train_pid_val in gpu_pids
            else:
                util = get_gpu_utilization()
                python_on_gpu = is_python_process_on_gpu()
                extend = (
                    (util is not None and util > gpu_busy_threshold) or python_on_gpu
                )
            if extend:
                deadline += extend_sec
                deadline = min(deadline, started + max_total_sec)
                if train_pid_val is not None:
                    print(
                        f"    [本任务训练进程 PID={train_pid_val} 仍在 GPU] 延长等待 {extend_sec}s，继续轮询…"
                    )
                else:
                    parts = []
                    if util is not None and util > gpu_busy_threshold:
                        parts.append(f"GPU 利用率 {util:.0f}%")
                    if python_on_gpu:
                        parts.append("检测到 Python 进程占用 GPU")
                    print(
                        f"    [{' | '.join(parts) or 'GPU 忙'}] 延长等待 {extend_sec}s，继续轮询…"
                    )
            else:
                reason = []
                if train_pid_val is not None:
                    reason.append(f"本任务训练进程 PID={train_pid_val} 已不在 GPU 上")
                else:
                    util = get_gpu_utilization()
                    if util is not None:
                        reason.append(f"GPU 利用率 {util}%")
                    reason.append("未检测到本任务或 Python 进程占用 GPU")
                raise RuntimeError(
                    f"Job {job_id} 未在 {max_wait_sec}s 内完成，当前 status={status}"
                    + (f"（{', '.join(reason)}）" if reason else "")
                )
        payload = request_json(
            base_url, "GET", f"/api/jobs/{job_id}", timeout_sec=timeout_sec
        )
        job = payload.get("job", {})
        status = str(job.get("status", ""))
        if status in ("success", "failed"):
            return job
        time.sleep(poll_interval_sec)


def download_asset(
    base_url: str, result_id: str, asset_name: str, out_path: Path
) -> None:
    url = f"{base_url.rstrip('/')}/api/results/{result_id}/assets/{asset_name}"
    req = Request(url=url, method="GET")
    with urlopen(req, timeout=60) as resp:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(resp.read())


def main() -> None:
    base_url = os.getenv("LABFLOW_BASE_URL", DEFAULT_BASE_URL)
    r_conda_env = os.getenv("LABFLOW_R_CONDA_ENV", DEFAULT_R_CONDA_ENV)
    r_conda_bat = os.getenv("LABFLOW_R_CONDA_BAT", DEFAULT_R_CONDA_BAT)
    # 默认用 TC.rds（筋膜）跑一条龙
    data_file = os.getenv("LABFLOW_FULLFLOW_DATA", "TC.rds")
    data_name = "TC-筋膜" if "TC.rds" in data_file else Path(data_file).stem

    rds_path = DATA_DIR / data_file
    if not rds_path.exists():
        print(f"数据文件不存在: {rds_path}")
        sys.exit(1)

    print("=== 全流程：转换 → 500×500 → 训练 → 预测 → 可视化报告 ===")
    print(f"BASE_URL={base_url}  DATA={data_file}")
    print(f"R: {r_conda_env}  conda.bat={r_conda_bat}")
    print()

    # ---------- 1. 注册 + 校验（R → h5ad）----------
    reg = request_json(
        base_url,
        "POST",
        "/api/datasets/register-local",
        payload={
            "local_path": str(rds_path),
            "dataset_name": data_name,
            "copy_to_upload_dir": True,
        },
        timeout_sec=60,
    )
    source_dataset_id = reg["dataset"]["id"]
    print(f"[1] 注册: dataset_id={source_dataset_id}")

    val = request_json(
        base_url,
        "POST",
        f"/api/datasets/{source_dataset_id}/validate",
        payload={
            "auto_convert": True,
            "r_exec_mode": "cmd_conda",
            "r_conda_env": r_conda_env,
            "r_conda_bat": r_conda_bat,
        },
        timeout_sec=600,
    )
    print(f"    校验/转换: status={val['dataset'].get('status')}")
    path_h5ad = val["dataset"].get("path_h5ad")
    if not path_h5ad or not Path(path_h5ad).exists():
        print("    path_h5ad 缺失，退出")
        sys.exit(1)

    # ---------- 2. Inspect + Prepare-training 500×500 ----------
    request_json(
        base_url,
        "POST",
        "/api/seurat/inspect",
        payload={"dataset_id": source_dataset_id, "umap_preview_limit": 500},
        timeout_sec=120,
    )
    group_col, cluster_col, selected = _infer_group_and_cluster_columns(Path(path_h5ad))
    prep = request_json(
        base_url,
        "POST",
        "/api/seurat/prepare-training",
        payload={
            "dataset_id": source_dataset_id,
            "group_column": group_col,
            "cluster_column": cluster_col,
            "selected_clusters": selected,
            "seed": 42,
        },
        timeout_sec=300,
    )
    prepared_dataset_id = prep["prepared_dataset_id"]
    n_cells = prep["n_cells"]
    n_genes = prep["n_genes"]
    print(
        f"[2] Prepare-training: prepared_dataset_id={prepared_dataset_id}  n_cells={n_cells}  n_genes={n_genes}"
    )

    # ---------- 3. 训练 ----------
    # 若后端为 LABFLOW_DRY_RUN=true，此处不会真实训练，仅写占位模型；真实模型需关闭 DRY_RUN 后重跑。
    train_payload = {
        "dataset_id": source_dataset_id,
        "prepared_dataset_id": prepared_dataset_id,
        "gene_size": n_genes,
        "output_dim": n_genes,
        "use_drug_structure": False,
        "batch_size": 32,
        "lr": 1e-4,
    }
    train_resp = request_json(
        base_url, "POST", "/api/jobs/train", payload=train_payload, timeout_sec=60
    )
    train_job_id = train_resp["job"]["id"]
    print(
        f"[3] 训练任务已提交: job_id={train_job_id}，轮询直至完成（超时后若 GPU 仍忙会自动延长）…"
    )
    train_job = poll_job_until_done(
        base_url,
        train_job_id,
        POLL_INTERVAL_SEC,
        TRAIN_POLL_TIMEOUT_SEC,
        extend_on_gpu_busy=True,
        gpu_busy_threshold=TRAIN_GPU_BUSY_THRESHOLD,
        extend_sec=TRAIN_EXTEND_SEC,
        max_total_sec=TRAIN_MAX_TOTAL_SEC,
    )
    if train_job.get("status") != "success":
        print(f"    训练失败: {train_job.get('error_msg')}")
        sys.exit(1)
    model_id = train_job.get("model_id")
    if not model_id:
        print("    训练成功但未返回 model_id")
        sys.exit(1)
    print(f"    训练完成: model_id={model_id}")
    art_dir = train_job.get("artifact_dir") or ""
    if art_dir:
        print(f"    模型/日志目录: {art_dir}（真实训练时 checkpoints/ 下为 .pt 模型）")

    # ---------- 4. 预测（用刚训练的模型 + prepared 数据）----------
    predict_payload = {
        "dataset_id": prepared_dataset_id,
        "model_id": model_id,
        "gene_size": n_genes,
        "output_dim": n_genes,
        "use_drug_structure": False,
    }
    predict_resp = request_json(
        base_url,
        "POST",
        "/api/jobs/predict",
        payload=predict_payload,
        timeout_sec=60,
    )
    predict_job_id = predict_resp["job"]["id"]
    print(f"[4] 预测任务已提交: job_id={predict_job_id}，轮询直至完成…")
    predict_job = poll_job_until_done(
        base_url, predict_job_id, POLL_INTERVAL_SEC, PREDICT_POLL_TIMEOUT_SEC
    )
    if predict_job.get("status") != "success":
        print(f"    预测失败: {predict_job.get('error_msg')}")
        sys.exit(1)
    result_id = predict_job.get("result_id")
    if not result_id:
        print("    预测成功但未返回 result_id")
        sys.exit(1)
    print(f"    预测完成: result_id={result_id}")

    # ---------- 5. 拉取结果与可视化报告 ----------
    result_resp = request_json(
        base_url,
        "GET",
        f"/api/results/job/{predict_job_id}",
        timeout_sec=30,
    )
    result = result_resp.get("result", {})
    summary = result.get("summary") or {}
    assets = summary.get("assets") or []
    print(
        f"[5] 可视化报告: n_cells={summary.get('n_cells')}  n_genes={summary.get('n_genes')}  assets={len(assets)}"
    )

    OUTPUT_REPORT_DIR.mkdir(parents=True, exist_ok=True)
    summary_out = OUTPUT_REPORT_DIR / "summary.json"
    summary_out.write_text(
        json.dumps(
            {
                "source_dataset_id": source_dataset_id,
                "prepared_dataset_id": prepared_dataset_id,
                "train_job_id": train_job_id,
                "model_id": model_id,
                "predict_job_id": predict_job_id,
                "result_id": result_id,
                "summary": summary,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(f"    已写入: {summary_out}")

    for i, item in enumerate(assets):
        if not isinstance(item, dict):
            continue
        path_str = item.get("path")
        name = (item.get("name") or Path(path_str).name) if path_str else f"asset_{i}"
        if not name:
            continue
        out_path = OUTPUT_REPORT_DIR / name
        try:
            download_asset(base_url, result_id, name, out_path)
            print(f"    已下载: {out_path}")
        except Exception as e:
            print(f"    下载 {name} 失败: {e}")

    print()
    print("=== 全流程结束 ===")
    print(f"报告目录: {OUTPUT_REPORT_DIR}")
    print("  summary.json  结果摘要")
    for f in OUTPUT_REPORT_DIR.iterdir():
        if f.suffix.lower() in (".png", ".jpg", ".jpeg"):
            print(f"  {f.name}  可视化图")


if __name__ == "__main__":
    main()
    sys.exit(0)
