"""LabFlow Windows launcher.

Purpose:
- Start backend (FastAPI/Uvicorn) and frontend (Vite) together on a Windows server.
- Support LAN access (bind to 0.0.0.0 by default).
- Work both as `python labflow_launcher.py` and PyInstaller-generated `.exe`.

This launcher assumes the project environment is already prepared:
- Python dependencies installed for backend
- Node dependencies installed in `frontend/`
"""

from __future__ import annotations

import argparse
import os
import shutil
import signal
import socket
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

WINDOWS = os.name == "nt"


@dataclass(frozen=True)
class LauncherConfig:
    project_root: Path
    uv_cmd: str
    npm_cmd: str
    host: str
    backend_port: int
    frontend_port: int
    frontend_mode: str  # preview | dev
    auto_build_frontend: bool
    backend_reload: bool


def _env_flag(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    value = raw.strip().lower()
    return value in {"1", "true", "yes", "on"}


def _ordered_unique(items: list[Path]) -> list[Path]:
    seen: set[str] = set()
    ordered: list[Path] = []
    for item in items:
        key = str(item.resolve())
        if key in seen:
            continue
        seen.add(key)
        ordered.append(item)
    return ordered


def _looks_like_project_root(path: Path) -> bool:
    return (path / "backend" / "app" / "main.py").exists() and (
        path / "frontend" / "package.json"
    ).exists()


def resolve_project_root(cli_project_root: str | None) -> Path:
    candidates: list[Path] = []

    if cli_project_root:
        candidates.append(Path(cli_project_root).expanduser())

    env_project_root = os.getenv("LABFLOW_PROJECT_ROOT")
    if env_project_root:
        candidates.append(Path(env_project_root).expanduser())

    if getattr(sys, "frozen", False):
        exe_dir = Path(sys.executable).resolve().parent
        candidates.extend([exe_dir, exe_dir.parent])

    script_dir = Path(__file__).resolve().parent
    candidates.extend([script_dir, Path.cwd()])

    expanded: list[Path] = []
    for base in _ordered_unique(candidates):
        expanded.append(base)
        expanded.extend(list(base.resolve().parents)[:4])

    for candidate in _ordered_unique(expanded):
        if _looks_like_project_root(candidate):
            return candidate

    checked = "\n  - " + "\n  - ".join(str(p) for p in _ordered_unique(expanded))
    raise RuntimeError(
        "Cannot locate LabFlow project root.\n"
        "Expected backend/app/main.py and frontend/package.json under one directory.\n"
        "Checked paths:" + checked + "\n"
        "Fix: place EXE in project root, or set LABFLOW_PROJECT_ROOT, or pass --project-root."
    )


def _command_exists(command: str) -> bool:
    if Path(command).exists():
        return True
    return shutil.which(command) is not None


def require_command(command: str, hint: str) -> None:
    if _command_exists(command):
        return
    raise RuntimeError(f"Missing command `{command}`. {hint}")


def detect_uv_command() -> str:
    env_cmd = os.getenv("LABFLOW_UV")
    if env_cmd:
        if not _command_exists(env_cmd):
            raise RuntimeError(
                f"LABFLOW_UV is set but command does not exist: {env_cmd}"
            )
        return env_cmd

    if _command_exists("uv"):
        return "uv"

    raise RuntimeError(
        "Missing command `uv`. Install uv, add it to PATH, or set LABFLOW_UV."
    )


def detect_npm_command() -> str:
    env_cmd = os.getenv("LABFLOW_NPM")
    if env_cmd:
        if not _command_exists(env_cmd):
            raise RuntimeError(
                f"LABFLOW_NPM is set but command does not exist: {env_cmd}"
            )
        return env_cmd

    npm_cmd = shutil.which("npm")
    if npm_cmd:
        return npm_cmd

    raise RuntimeError(
        "Missing command `npm`. Install Node.js, add npm to PATH, or set LABFLOW_NPM."
    )


def get_local_ip() -> str:
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        sock.connect(("8.8.8.8", 80))
        return sock.getsockname()[0]
    except Exception:
        return "127.0.0.1"
    finally:
        sock.close()


def wait_http_ready(url: str, timeout_seconds: int) -> bool:
    deadline = time.time() + timeout_seconds
    while time.time() < deadline:
        try:
            with urlopen(url, timeout=2) as response:  # noqa: S310
                if 200 <= response.status < 500:
                    return True
        except URLError:
            pass
        except Exception:
            pass
        time.sleep(1)
    return False


def run_checked(cmd: list[str], cwd: Path, name: str) -> None:
    print(f"[run] {name}: {' '.join(cmd)}")
    proc = subprocess.run(cmd, cwd=str(cwd), check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"{name} failed with exit code {proc.returncode}")


def spawn_process(
    cmd: list[str],
    cwd: Path,
    name: str,
    env: dict[str, str] | None = None,
) -> subprocess.Popen:
    print(f"[start] {name}: {' '.join(cmd)}")
    kwargs: dict[str, object] = {
        "cwd": str(cwd),
        "env": env,
    }

    if WINDOWS:
        kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP  # type: ignore[attr-defined]
    else:
        kwargs["start_new_session"] = True

    return subprocess.Popen(cmd, **kwargs)  # noqa: S603


def stop_process(proc: subprocess.Popen | None, name: str) -> None:
    if proc is None or proc.poll() is not None:
        return

    print(f"[stop] {name} (pid={proc.pid})")
    try:
        if WINDOWS:
            proc.send_signal(signal.CTRL_BREAK_EVENT)  # type: ignore[attr-defined]
            proc.wait(timeout=8)
        else:
            proc.terminate()
            proc.wait(timeout=8)
        return
    except Exception:
        pass

    print(f"[stop] force-kill {name} (pid={proc.pid})")
    if WINDOWS:
        subprocess.run(
            ["taskkill", "/PID", str(proc.pid), "/T", "/F"],
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    else:
        proc.kill()


def ensure_project_layout(project_root: Path) -> None:
    required_paths = [
        project_root / "backend" / "app" / "main.py",
        project_root / "frontend" / "package.json",
    ]
    for path in required_paths:
        if not path.exists():
            raise RuntimeError(f"Required file not found: {path}")


def prepare_frontend_if_needed(config: LauncherConfig) -> None:
    frontend_dir = config.project_root / "frontend"
    if not (frontend_dir / "node_modules").exists():
        raise RuntimeError(
            "frontend/node_modules not found. Run `cd frontend && npm install` first."
        )

    if config.frontend_mode != "preview":
        return

    dist_dir = frontend_dir / "dist"
    if config.auto_build_frontend or (not dist_dir.exists()):
        run_checked([config.npm_cmd, "run", "build"], frontend_dir, "frontend build")


def build_backend_cmd(config: LauncherConfig) -> list[str]:
    cmd = [
        config.uv_cmd,
        "run",
        "--active",
        "--no-project",
        "python",
        "-m",
        "uvicorn",
        "backend.app.main:app",
        "--host",
        config.host,
        "--port",
        str(config.backend_port),
    ]
    if config.backend_reload:
        cmd.append("--reload")
    return cmd


def build_frontend_cmd(config: LauncherConfig) -> list[str]:
    base = [
        config.npm_cmd,
        "run",
        "preview" if config.frontend_mode == "preview" else "dev",
    ]
    return [
        *base,
        "--",
        "--host",
        config.host,
        "--port",
        str(config.frontend_port),
        "--strictPort",
    ]


def print_banner(config: LauncherConfig) -> None:
    local_ip = get_local_ip()
    print("=" * 76)
    print("LabFlow Windows Launcher")
    print(f"project root: {config.project_root}")
    print(f"frontend mode: {config.frontend_mode}")
    print(f"backend:  http://localhost:{config.backend_port}")
    print(f"backend (LAN):  http://{local_ip}:{config.backend_port}")
    print(f"frontend: http://localhost:{config.frontend_port}")
    print(f"frontend (LAN): http://{local_ip}:{config.frontend_port}")
    print("=" * 76)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Start LabFlow backend+frontend on Windows for LAN usage."
    )
    parser.add_argument(
        "--project-root", default=None, help="Path to LabFlow project root"
    )
    parser.add_argument(
        "--host",
        default=os.getenv("LABFLOW_HOST", "0.0.0.0"),
        help="Bind host for backend/frontend (default: 0.0.0.0)",
    )
    parser.add_argument(
        "--backend-port",
        type=int,
        default=int(os.getenv("LABFLOW_BACKEND_PORT", "8000")),
        help="Backend port (default: 8000)",
    )
    parser.add_argument(
        "--frontend-port",
        type=int,
        default=int(os.getenv("LABFLOW_FRONTEND_PORT", "5173")),
        help="Frontend port (default: 5173)",
    )
    parser.add_argument(
        "--frontend-mode",
        choices=["preview", "dev"],
        default=os.getenv("LABFLOW_FRONTEND_MODE", "preview"),
        help="Use preview for server-like runtime (default) or dev for debugging.",
    )
    parser.add_argument(
        "--build-frontend",
        action="store_true",
        help="Force run `npm run build` before starting frontend.",
    )
    parser.add_argument(
        "--backend-reload",
        action="store_true",
        help="Enable uvicorn --reload (for development only).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved config and commands without starting processes.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        project_root = resolve_project_root(args.project_root)
        ensure_project_layout(project_root)

        uv_cmd = detect_uv_command()
        npm_cmd = detect_npm_command()

        config = LauncherConfig(
            project_root=project_root,
            uv_cmd=uv_cmd,
            npm_cmd=npm_cmd,
            host=args.host,
            backend_port=args.backend_port,
            frontend_port=args.frontend_port,
            frontend_mode=args.frontend_mode,
            auto_build_frontend=args.build_frontend
            or _env_flag("LABFLOW_AUTO_BUILD_FRONTEND", False),
            backend_reload=args.backend_reload
            or _env_flag("LABFLOW_BACKEND_RELOAD", False),
        )

        print_banner(config)

        backend_cmd = build_backend_cmd(config)
        frontend_cmd = build_frontend_cmd(config)

        if args.dry_run:
            print("[dry-run] backend command:")
            print("  " + " ".join(backend_cmd))
            print("[dry-run] frontend command:")
            print("  " + " ".join(frontend_cmd))
            return 0

        prepare_frontend_if_needed(config)

        backend_proc: subprocess.Popen | None = None
        frontend_proc: subprocess.Popen | None = None

        try:
            backend_proc = spawn_process(
                backend_cmd,
                config.project_root,
                "backend",
            )

            backend_health = f"http://127.0.0.1:{config.backend_port}/api/health"
            if not wait_http_ready(backend_health, timeout_seconds=30):
                raise RuntimeError(
                    "Backend did not become healthy within 30s: " + backend_health
                )
            print("[ok] backend health check passed")

            frontend_env = os.environ.copy()
            vite_api_base = os.getenv("LABFLOW_VITE_API_BASE")
            if vite_api_base:
                frontend_env["VITE_API_BASE"] = vite_api_base

            frontend_proc = spawn_process(
                frontend_cmd,
                config.project_root / "frontend",
                "frontend",
                env=frontend_env,
            )

            frontend_url = f"http://127.0.0.1:{config.frontend_port}"
            if wait_http_ready(frontend_url, timeout_seconds=45):
                print("[ok] frontend is reachable")
            else:
                print("[warn] frontend not reachable yet; it may still be compiling")

            print("[info] Press Ctrl+C in this window to stop both services.")
            while True:
                time.sleep(1)
                if backend_proc.poll() is not None:
                    raise RuntimeError("Backend process exited unexpectedly.")
                if frontend_proc.poll() is not None:
                    raise RuntimeError("Frontend process exited unexpectedly.")

        except KeyboardInterrupt:
            print("\n[info] Received Ctrl+C, shutting down...")
        finally:
            stop_process(frontend_proc, "frontend")
            stop_process(backend_proc, "backend")
            print("[done] LabFlow services stopped.")

    except Exception as exc:  # noqa: BLE001
        print(f"[error] {exc}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
