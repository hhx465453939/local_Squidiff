"""
LabFlow Windows 启动器（实验室内网一键启动）

功能：
- 一键启动 LabFlow 后端（FastAPI + Uvicorn）与前端（React + Vite）
- 默认监听 0.0.0.0，局域网其他机器可通过本机 IP 访问
- 控制台中展示本机访问地址与内网访问地址

使用前提（必须先完成，否则启动会失败）：
- 已按 README「5. 开发与部署」准备好 Python 环境，并安装好后端依赖
- 已在 frontend 目录执行过一次 `npm install`
- 在 Windows 上运行（PowerShell / CMD 双击均可）

打包为 .exe 的方式请见 `docs/Windows一键启动器.md`
"""

from __future__ import annotations

import os
import signal
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional


PROJECT_ROOT = Path(__file__).resolve().parent
FRONTEND_DIR = PROJECT_ROOT / "frontend"


def detect_python_command() -> str:
    """
    检测用于启动后端的 Python 命令。

    优先级：
    1. 环境变量 LABFLOW_PYTHON
    2. 'python'
    """
    env_cmd = os.environ.get("LABFLOW_PYTHON")
    if env_cmd:
        return env_cmd
    return "python"


def get_local_ip() -> str:
    """
    获取当前主机在局域网中的 IPv4 地址。

    若无法自动检测，则退回到 127.0.0.1。
    """
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        try:
            # 这个地址不会真的被访问，只是用于触发本机选路
            s.connect(("8.8.8.8", 80))
            ip = s.getsockname()[0]
        finally:
            s.close()
        return ip
    except Exception:
        return "127.0.0.1"


def spawn_process(
    cmd: List[str],
    cwd: Optional[Path] = None,
    env: Optional[dict] = None,
    name: str = "",
) -> subprocess.Popen:
    """
    启动子进程（不截获输出，直接继承当前控制台）。
    """
    display_name = name or " ".join(cmd)
    print(f"[启动] {display_name}")
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd else None,
            env=env,
        )
    except FileNotFoundError:
        print(f"[错误] 无法找到命令：{cmd[0]}，请确认已安装并加入 PATH。")
        raise
    return proc


def start_backend(python_cmd: str) -> subprocess.Popen:
    """
    启动 Uvicorn 后端（监听 0.0.0.0:8000）。
    """
    cmd = [
        python_cmd,
        "-m",
        "uvicorn",
        "backend.app.main:app",
        "--host",
        "0.0.0.0",
        "--port",
        "8000",
    ]
    return spawn_process(cmd, cwd=PROJECT_ROOT, name="LabFlow 后端 (Uvicorn)")


def start_frontend() -> subprocess.Popen:
    """
    启动 Vite 前端开发服务器（监听 0.0.0.0:5173）。

    Vite 的 host/port 已在 frontend/vite.config.ts 中配置。
    """
    cmd = ["npm", "run", "dev"]
    env = os.environ.copy()
    # 如需指定后端地址，可在此处设置 VITE_API_BASE
    return spawn_process(cmd, cwd=FRONTEND_DIR, env=env, name="LabFlow 前端 (Vite)")


def print_banner():
    print("=" * 70)
    print("LabFlow Windows 启动器")
    print("- 一键启动后端（FastAPI + Uvicorn）和前端（React + Vite）")
    print("- 默认监听 0.0.0.0，支持本机与内网访问")
    print("=" * 70)
    print()


def print_urls():
    ip = get_local_ip()
    print("访问方式：")
    print(f"- 本机浏览：  http://localhost:5173")
    print(f"- 内网访问：  http://{ip}:5173")
    print()
    print("后端 API：")
    print(f"- 本机浏览：  http://localhost:8000")
    print(f"- 内网访问：  http://{ip}:8000")
    print()
    print("提示：")
    print("- 首次使用请先按照 README 完成 Python/Node 环境与依赖安装；")
    print("- 首次在 frontend 目录运行 `npm install` 可能耗时较长；")
    print("- 若看到 Uvicorn 和 Vite 启动日志，即表示服务已启动；")
    print("- 若要停止服务，请在此窗口按 Ctrl+C。")
    print()


def ensure_project_structure():
    """
    简单检查项目结构是否完整。
    """
    backend_main = PROJECT_ROOT / "backend" / "app" / "main.py"
    if not backend_main.exists():
        print("[错误] 未找到 backend/app/main.py，请确认在 LabFlow 项目根目录运行本程序。")
        sys.exit(1)

    if not FRONTEND_DIR.exists():
        print("[错误] 未找到 frontend 目录，请确认在 LabFlow 项目根目录运行本程序。")
        sys.exit(1)


def main():
    print_banner()
    ensure_project_structure()

    python_cmd = detect_python_command()
    print(f"[信息] 使用 Python 命令：{python_cmd}")
    print()

    backend_proc: Optional[subprocess.Popen] = None
    frontend_proc: Optional[subprocess.Popen] = None

    try:
        backend_proc = start_backend(python_cmd)
        # 稍等片刻再启动前端，避免日志混在一起难以分辨
        time.sleep(2)
        frontend_proc = start_frontend()

        print()
        print_urls()

        # 主循环：仅用于保持进程存活，并响应 Ctrl+C
        while True:
            time.sleep(1)
            # 若任一子进程退出，则提示用户
            if backend_proc.poll() is not None:
                print("[警告] 后端进程已退出，请检查上方日志。")
                break
            if frontend_proc.poll() is not None:
                print("[警告] 前端进程已退出，请检查上方日志。")
                break

    except KeyboardInterrupt:
        print("\n[信息] 收到中断信号，正在关闭子进程...")
    finally:
        # 依次尝试优雅关闭子进程
        for name, proc in [
            ("前端 (Vite)", frontend_proc),
            ("后端 (Uvicorn)", backend_proc),
        ]:
            if proc and proc.poll() is None:
                print(f"[停止] 终止 {name} 进程...")
                try:
                    if os.name == "nt":
                        proc.send_signal(signal.CTRL_BREAK_EVENT)  # type: ignore[attr-defined]
                    else:
                        proc.terminate()
                except Exception:
                    proc.terminate()

                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    print(f"[强制] {name} 未在 10 秒内退出，尝试强制杀死...")
                    proc.kill()

        print("[完成] LabFlow 服务已停止，可以关闭窗口。")


if __name__ == "__main__":
    main()

