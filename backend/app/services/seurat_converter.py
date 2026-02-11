from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from ..core.config import settings


SUPPORTED_INPUTS = {".h5ad", ".h5seurat", ".rds"}


def _build_r_command(
    script_path: Path,
    input_path: Path,
    output_path: Path,
    *,
    r_exec_mode: str | None,
    r_conda_env: str | None,
    r_conda_bat: str | None,
    rscript_bin: str | None,
) -> tuple[list[str], Path | None]:
    """Build R invocation command with optional Windows cmd+conda activation.
    Returns (argv, tmp_bat_path). If tmp_bat_path is set, caller must unlink it after run.
    """
    exec_mode = r_exec_mode or settings.r_exec_mode
    conda_env = r_conda_env or settings.r_conda_env
    conda_bat = r_conda_bat or settings.r_conda_bat or "conda.bat"
    rscript = rscript_bin or settings.rscript_bin

    if exec_mode == "cmd_conda":
        if not conda_env:
            raise RuntimeError(
                "LABFLOW_R_EXEC_MODE=cmd_conda requires LABFLOW_R_CONDA_ENV"
            )
        # Windows: pass conda activation + R via a temp .bat to avoid cmd /c quoting issues.
        # Use forward slashes for R paths so backslashes are not interpreted as escapes in R.
        script_s = str(script_path).replace("\\", "/")
        input_s = str(input_path).replace("\\", "/")
        output_s = str(output_path).replace("\\", "/")
        tmp = tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".bat",
            delete=False,
            encoding="utf-8",
        )
        tmp.write(
            f"@echo off\n"
            f'call "{conda_bat}" activate "{conda_env}"\n'
            f'"{rscript}" "{script_s}" "{input_s}" "{output_s}"\n'
        )
        tmp.close()
        return ["cmd", "/c", tmp.name], Path(tmp.name)

    if shutil.which(rscript) is None and not Path(rscript).exists():
        raise RuntimeError(
            "Rscript was not found. Please install R and ensure Rscript is in PATH, "
            "or set LABFLOW_RSCRIPT_BIN / rscript_bin to the full path."
        )

    return [
        rscript,
        str(script_path),
        str(input_path),
        str(output_path),
    ], None


def convert_to_h5ad(
    input_path: Path,
    output_dir: Path,
    r_exec_mode: str | None = None,
    r_conda_env: str | None = None,
    r_conda_bat: str | None = None,
    rscript_bin: str | None = None,
) -> Path:
    """Convert Seurat input to h5ad via an R helper script."""
    suffix = input_path.suffix.lower()
    if suffix not in SUPPORTED_INPUTS:
        raise ValueError(f"Unsupported input format: {suffix}")

    if suffix == ".h5ad":
        return input_path

    script_path = Path(__file__).resolve().parents[2] / "scripts" / "seurat_to_h5ad.R"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{input_path.stem}.h5ad"

    cmd, tmp_bat = _build_r_command(
        script_path,
        input_path,
        output_path,
        r_exec_mode=r_exec_mode,
        r_conda_env=r_conda_env,
        r_conda_bat=r_conda_bat,
        rscript_bin=rscript_bin,
    )
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                "Seurat conversion failed. "
                f"cmd={' '.join(cmd)} "
                f"stdout={proc.stdout.strip()} stderr={proc.stderr.strip()}"
            )
    finally:
        if tmp_bat is not None and tmp_bat.exists():
            try:
                tmp_bat.unlink()
            except OSError:
                pass

    if not output_path.exists():
        extra = ""
        if proc.stdout.strip() or proc.stderr.strip():
            extra = f" stdout={proc.stdout.strip()!r} stderr={proc.stderr.strip()!r}"
        raise RuntimeError(
            f"Expected converted file not found: {output_path}. R exited 0.{extra}"
        )
    return output_path
