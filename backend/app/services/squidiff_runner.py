from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np

from ..core.config import settings
from .visualize import generate_visual_assets


class SquidiffRunner:
    def __init__(self, repo_root: Path) -> None:
        self.repo_root = repo_root

    def _write_log(self, log_path: Path, content: str) -> None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(content)
            if not content.endswith("\n"):
                handle.write("\n")

    def run_train(
        self,
        *,
        job_dir: Path,
        params: dict[str, Any],
        on_start: Callable[[int], None] | None = None,
    ) -> dict[str, Any]:
        job_dir.mkdir(parents=True, exist_ok=True)
        log_path = job_dir / "train.log"
        checkpoint_dir = job_dir / "checkpoints"
        checkpoint_dir.mkdir(parents=True, exist_ok=True)

        if settings.dry_run:
            fake_model = checkpoint_dir / "model.pt"
            fake_model.write_bytes(b"dry-run-model")
            self._write_log(log_path, "[dry-run] train completed")
            return {
                "model_path": str(fake_model),
                "log_path": str(log_path),
                "metrics": {"dry_run": True},
            }

        cmd = [
            sys.executable,
            str(self.repo_root / "train_squidiff.py"),
            "--logger_path",
            str(job_dir / "logger"),
            "--data_path",
            params["data_path"],
            "--resume_checkpoint",
            str(checkpoint_dir),
            "--gene_size",
            str(params["gene_size"]),
            "--output_dim",
            str(params["output_dim"]),
            "--batch_size",
            str(params.get("batch_size", 64)),
            "--lr",
            str(params.get("lr", 1e-4)),
        ]
        if params.get("use_drug_structure"):
            cmd.extend(["--use_drug_structure", "True"])
            if params.get("control_data_path"):
                cmd.extend(["--control_data_path", params["control_data_path"]])

        proc = subprocess.Popen(
            cmd,
            cwd=str(self.repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if on_start is not None:
            on_start(proc.pid)
        stdout, stderr = proc.communicate()
        self._write_log(log_path, stdout or "")
        self._write_log(log_path, stderr or "")

        if proc.returncode != 0:
            tail = (stderr or "").strip() or (stdout or "").strip()
            if not tail and log_path.exists():
                try:
                    tail = log_path.read_text(encoding="utf-8", errors="replace").strip()
                except OSError:
                    pass
            if len(tail) > 1500:
                tail = tail[-1500:]
            raise RuntimeError(
                f"Training command failed with exit code {proc.returncode}"
                + (f". Last output:\n{tail}" if tail else "")
            )

        model_path = self._discover_model_path(checkpoint_dir)
        return {
            "model_path": str(model_path),
            "log_path": str(log_path),
            "metrics": {"dry_run": False},
        }

    def run_predict(self, *, job_dir: Path, params: dict[str, Any]) -> dict[str, Any]:
        job_dir.mkdir(parents=True, exist_ok=True)
        log_path = job_dir / "predict.log"
        result_dir = job_dir / "result"
        result_dir.mkdir(parents=True, exist_ok=True)

        if settings.dry_run:
            fake = np.random.default_rng(42).normal(size=(128, params["gene_size"]))
            pred_path = result_dir / "prediction.npy"
            np.save(pred_path, fake)
            summary = generate_visual_assets(
                prediction_matrix=fake,
                output_dir=result_dir / "viz",
            )
            self._write_log(log_path, "[dry-run] prediction completed")
            return {
                "prediction_path": str(pred_path),
                "log_path": str(log_path),
                "summary": summary,
            }

        data_path = Path(params["data_path"])
        model_path = Path(params["model_path"])
        import scanpy as sc

        adata = sc.read_h5ad(str(data_path))
        matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
        matrix = np.asarray(matrix, dtype=np.float32)

        # Import only at runtime to keep startup light.
        import torch
        import sample_squidiff

        device = "cuda" if torch.cuda.is_available() else "cpu"
        sampler = sample_squidiff.sampler(
            model_path=str(model_path),
            gene_size=int(params["gene_size"]),
            output_dim=int(params["output_dim"]),
            use_drug_structure=bool(params.get("use_drug_structure", False)),
        )

        tensor_x = torch.tensor(matrix, dtype=torch.float32).to(device)
        with torch.no_grad():
            z_sem = sampler.model.encoder(tensor_x)
            pred = sampler.pred(z_sem, gene_size=matrix.shape[1])
        pred_np = pred.detach().cpu().numpy()

        pred_path = result_dir / "prediction.npy"
        np.save(pred_path, pred_np)

        summary = generate_visual_assets(
            prediction_matrix=pred_np,
            output_dir=result_dir / "viz",
        )
        self._write_log(
            log_path,
            json.dumps(
                {"prediction_path": str(pred_path), "summary": summary}, indent=2
            ),
        )
        return {
            "prediction_path": str(pred_path),
            "log_path": str(log_path),
            "summary": summary,
        }

    def _discover_model_path(self, checkpoint_dir: Path) -> Path:
        candidates = sorted(checkpoint_dir.rglob("*.pt"))
        if candidates:
            return candidates[-1]
        fallback = checkpoint_dir / "model.pt"
        if not fallback.exists():
            raise RuntimeError(
                "Training completed but no model checkpoint was found under "
                f"{checkpoint_dir}"
            )
        return fallback
