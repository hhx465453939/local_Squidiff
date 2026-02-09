from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


@dataclass
class DatasetUatResult:
    dataset_id: str
    passed: bool
    inspect_ok: bool
    prepare_ok: bool
    train_ok: bool
    prepared_dataset_id: str | None
    n_cells: int | None
    n_genes: int | None
    train_job_id: str | None
    train_status: str | None
    detail: str


def request_json(
    *,
    base_url: str,
    method: str,
    path: str,
    payload: dict[str, Any] | None = None,
    timeout_sec: int = 30,
) -> dict[str, Any]:
    url = f"{base_url.rstrip('/')}{path}"
    headers = {"Content-Type": "application/json"}
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    request = Request(url=url, data=data, headers=headers, method=method.upper())
    try:
        with urlopen(request, timeout=timeout_sec) as response:
            body = response.read().decode("utf-8")
    except HTTPError as exc:
        error_body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(
            f"{method.upper()} {path} failed: {exc.code} {error_body}"
        ) from exc
    except URLError as exc:
        raise RuntimeError(f"{method.upper()} {path} failed: {exc.reason}") from exc

    try:
        return json.loads(body)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            f"{method.upper()} {path} returned non-JSON: {body}"
        ) from exc


def poll_train_job_until_done(
    *,
    base_url: str,
    job_id: str,
    poll_interval_sec: int,
    max_wait_sec: int,
    timeout_sec: int,
) -> dict[str, Any]:
    deadline = time.time() + max_wait_sec
    while time.time() <= deadline:
        payload = request_json(
            base_url=base_url,
            method="GET",
            path=f"/api/jobs/{job_id}",
            timeout_sec=timeout_sec,
        )
        job = payload.get("job", {})
        status = str(job.get("status", ""))
        if status in {"success", "failed"}:
            return job
        time.sleep(poll_interval_sec)
    raise RuntimeError(f"train job timeout after {max_wait_sec}s: {job_id}")


def run_dataset_uat(
    *,
    base_url: str,
    dataset_id: str,
    group_column: str,
    cluster_column: str,
    selected_clusters: list[str],
    seed: int,
    skip_train: bool,
    poll_interval_sec: int,
    max_wait_sec: int,
    timeout_sec: int,
) -> DatasetUatResult:
    inspect_ok = False
    prepare_ok = False
    train_ok = False
    prepared_dataset_id: str | None = None
    n_cells: int | None = None
    n_genes: int | None = None
    train_job_id: str | None = None
    train_status: str | None = None

    try:
        inspect_payload = request_json(
            base_url=base_url,
            method="POST",
            path="/api/seurat/inspect",
            payload={"dataset_id": dataset_id, "umap_preview_limit": 300},
            timeout_sec=timeout_sec,
        )
        inspect = inspect_payload.get("inspect")
        if not isinstance(inspect, dict):
            raise RuntimeError("inspect response missing 'inspect' object")
        metadata_columns = inspect.get("metadata_columns")
        if not isinstance(metadata_columns, list):
            raise RuntimeError("inspect metadata_columns is invalid")
        inspect_ok = True

        prepare_payload = request_json(
            base_url=base_url,
            method="POST",
            path="/api/seurat/prepare-training",
            payload={
                "dataset_id": dataset_id,
                "group_column": group_column,
                "cluster_column": cluster_column,
                "selected_clusters": selected_clusters,
                "seed": seed,
            },
            timeout_sec=timeout_sec,
        )
        prepared_dataset_id = str(prepare_payload.get("prepared_dataset_id", ""))
        n_cells = int(prepare_payload.get("n_cells", -1))
        n_genes = int(prepare_payload.get("n_genes", -1))

        if not prepared_dataset_id:
            raise RuntimeError("prepare-training response missing prepared_dataset_id")
        if n_cells <= 0 or n_cells > 500:
            raise RuntimeError(f"n_cells out of range: {n_cells}")
        if n_genes <= 0 or n_genes > 500:
            raise RuntimeError(f"n_genes out of range: {n_genes}")
        prepare_ok = True

        if skip_train:
            train_ok = True
            train_status = "skipped"
        else:
            train_submit = request_json(
                base_url=base_url,
                method="POST",
                path="/api/jobs/train",
                payload={
                    "dataset_id": dataset_id,
                    "prepared_dataset_id": prepared_dataset_id,
                    "gene_size": n_genes,
                    "output_dim": n_genes,
                    "use_drug_structure": False,
                    "batch_size": 32,
                    "lr": 1e-4,
                },
                timeout_sec=timeout_sec,
            )
            job = train_submit.get("job")
            if not isinstance(job, dict) or not isinstance(job.get("id"), str):
                raise RuntimeError("train submit response missing job.id")
            train_job_id = job["id"]

            final_job = poll_train_job_until_done(
                base_url=base_url,
                job_id=train_job_id,
                poll_interval_sec=poll_interval_sec,
                max_wait_sec=max_wait_sec,
                timeout_sec=timeout_sec,
            )
            train_status = str(final_job.get("status"))
            train_ok = train_status == "success"
            if not train_ok:
                error_text = final_job.get("error_msg", "")
                raise RuntimeError(f"train job failed: {train_job_id} {error_text}")

        return DatasetUatResult(
            dataset_id=dataset_id,
            passed=inspect_ok and prepare_ok and train_ok,
            inspect_ok=inspect_ok,
            prepare_ok=prepare_ok,
            train_ok=train_ok,
            prepared_dataset_id=prepared_dataset_id,
            n_cells=n_cells,
            n_genes=n_genes,
            train_job_id=train_job_id,
            train_status=train_status,
            detail="ok",
        )
    except Exception as exc:  # noqa: BLE001
        return DatasetUatResult(
            dataset_id=dataset_id,
            passed=False,
            inspect_ok=inspect_ok,
            prepare_ok=prepare_ok,
            train_ok=train_ok,
            prepared_dataset_id=prepared_dataset_id,
            n_cells=n_cells,
            n_genes=n_genes,
            train_job_id=train_job_id,
            train_status=train_status,
            detail=str(exc),
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Phase 4 UAT runner for Seurat V2 pipeline"
    )
    parser.add_argument(
        "--base-url",
        default="http://localhost:8000",
        help="Backend base URL",
    )
    parser.add_argument(
        "--dataset-id",
        action="append",
        dest="dataset_ids",
        required=True,
        help="Dataset ID to test (repeat at least twice)",
    )
    parser.add_argument("--group-column", required=True, help="group_column value")
    parser.add_argument("--cluster-column", required=True, help="cluster_column value")
    parser.add_argument(
        "--selected-clusters",
        required=True,
        help="Comma separated selected cluster values, e.g. T,B,NK",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--skip-train",
        action="store_true",
        help="Skip train job execution and only verify inspect + prepare",
    )
    parser.add_argument(
        "--poll-interval-sec",
        type=int,
        default=5,
        help="Polling interval when waiting for train job",
    )
    parser.add_argument(
        "--max-wait-sec",
        type=int,
        default=1800,
        help="Max wait seconds for one train job",
    )
    parser.add_argument("--timeout-sec", type=int, default=30, help="HTTP timeout")
    parser.add_argument(
        "--report-path",
        default="tmp_cache/uat_phase4_seurat_v2_report.json",
        help="Output report path",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    dataset_ids: list[str] = args.dataset_ids
    if len(dataset_ids) < 2:
        print("ERROR: at least two --dataset-id values are required for Phase 4 UAT.")
        return 2

    selected_clusters = [
        item.strip() for item in str(args.selected_clusters).split(",") if item.strip()
    ]
    if not selected_clusters:
        print("ERROR: --selected-clusters must include at least one non-empty value.")
        return 2

    results: list[DatasetUatResult] = []
    for dataset_id in dataset_ids:
        result = run_dataset_uat(
            base_url=args.base_url,
            dataset_id=dataset_id,
            group_column=args.group_column,
            cluster_column=args.cluster_column,
            selected_clusters=selected_clusters,
            seed=args.seed,
            skip_train=bool(args.skip_train),
            poll_interval_sec=args.poll_interval_sec,
            max_wait_sec=args.max_wait_sec,
            timeout_sec=args.timeout_sec,
        )
        results.append(result)
        status = "PASS" if result.passed else "FAIL"
        print(
            f"[{status}] dataset={result.dataset_id} "
            f"inspect={result.inspect_ok} prepare={result.prepare_ok} train={result.train_ok} "
            f"prepared_dataset_id={result.prepared_dataset_id or '-'} detail={result.detail}"
        )

    passed_count = sum(1 for item in results if item.passed)
    report = {
        "timestamp": int(time.time()),
        "base_url": args.base_url,
        "dataset_ids": dataset_ids,
        "group_column": args.group_column,
        "cluster_column": args.cluster_column,
        "selected_clusters": selected_clusters,
        "seed": args.seed,
        "skip_train": bool(args.skip_train),
        "passed_count": passed_count,
        "total_count": len(results),
        "results": [asdict(item) for item in results],
    }
    report_path = Path(args.report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"Report saved: {report_path}")

    if passed_count != len(results):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
