# Windows 一键启动器（内网部署版）

目标：在 Windows 服务器上“一次部署后双击启动”，同时拉起后端与前端，供实验室内网用户访问。

---

## 1. 适用场景

- 一台固定 Windows 机器作为实验室内网服务机。
- LabFlow 代码和依赖已提前部署完成。
- 需要非技术同学通过浏览器直接使用系统。

---

## 2. 启动器能力（当前版本）

`labflow_launcher.py` / `LabFlowLauncher.exe` 会：

1. 自动定位项目根目录（支持脚本运行和 EXE 运行）。
2. 启动后端：`uv run --active --no-project python -m uvicorn backend.app.main:app --host 0.0.0.0 --port 8000`。
3. 启动前端：
- 默认 `preview` 模式（适合服务器长期运行）。
- 启动命令：`npm run preview -- --host 0.0.0.0 --port 5173 --strictPort`。
4. 打印本机与内网访问地址。
5. Ctrl+C 时同时关闭前后端进程。

---

## 3. 前置条件（只需做一次）

1. 完成项目部署（见 `README.md` 与 `docs/部署文档.md`）。
2. 已安装并可用：
- Python（建议 3.11）
- Node.js / npm
3. 已安装后端依赖与前端依赖：
```powershell
pip install -r requirements.txt -r backend/requirements.txt
cd frontend
npm install
```

---

## 4. 开发者本地运行（脚本）

在项目根目录执行：

```powershell
python labflow_launcher.py
```

如果只想检查配置不真正启动：

```powershell
python labflow_launcher.py --dry-run
```

---

## 5. 打包 EXE（给运维/管理员）

```powershell
pip install pyinstaller
python -m PyInstaller --onefile --name LabFlowLauncher labflow_launcher.py
```

输出文件：
- `dist/LabFlowLauncher.exe`

推荐：
- 将 EXE 放在项目根目录同级（最省事）。
- 或设置 `LABFLOW_PROJECT_ROOT` 指向项目根目录。

---

## 6. EXE 启动与内网访问

双击 `LabFlowLauncher.exe` 后，服务默认监听：

- 前端：`0.0.0.0:5173`
- 后端：`0.0.0.0:8000`

访问方式：
- 本机：`http://localhost:5173`
- 内网用户：`http://<服务器IP>:5173`

注意：
- 防火墙需放行 5173 与 8000。
- 若前端页面能开但 API 报错，先检查后端 `8000` 是否启动成功。

---

## 7. 常用参数与环境变量

### 7.1 命令行参数

```powershell
python labflow_launcher.py `
  --project-root E:\Development\local_Squidiff `
  --host 0.0.0.0 `
  --backend-port 8000 `
  --frontend-port 5173 `
  --frontend-mode preview
```

可选：
- `--frontend-mode preview|dev`（默认 `preview`）
- `--build-frontend`（强制先执行 `npm run build`）
- `--backend-reload`（开发调试用）

### 7.2 环境变量

- `LABFLOW_PROJECT_ROOT`：项目根目录
- `LABFLOW_UV`：指定 `uv` 可执行文件（当 `uv` 不在 PATH 时）
- `LABFLOW_NPM`：指定 `npm` 可执行文件（当 `npm` 不在 PATH 或启动时报 WinError 2 时）
- `LABFLOW_HOST`：监听地址（默认 `0.0.0.0`）
- `LABFLOW_BACKEND_PORT`：后端端口（默认 `8000`）
- `LABFLOW_FRONTEND_PORT`：前端端口（默认 `5173`）
- `LABFLOW_FRONTEND_MODE`：`preview` 或 `dev`
- `LABFLOW_AUTO_BUILD_FRONTEND`：`true/false`
- `LABFLOW_BACKEND_RELOAD`：`true/false`
- `LABFLOW_VITE_API_BASE`：覆盖前端 API 地址（一般不需要）

---

## 8. 常见问题

1. 启动器提示找不到 Python
- 设置 `LABFLOW_PYTHON` 为完整路径（例如 `F:\software\Miniconda3\python.exe`）。

2. 启动器提示找不到 npm
- 安装 Node.js 并确认 `npm -v` 可用。

3. 启动器提示找不到 `uv`
- 安装 uv 并确认 `uv --version` 可用，或设置 `LABFLOW_UV` 为完整路径（例如 `C:\Users\<用户名>\AppData\Roaming\Python\Python312\Scripts\uv.exe`）。

4. 启动前端时报 `[WinError 2] 系统找不到指定的文件`
- 通常是 `npm` 命令解析失败。设置 `LABFLOW_NPM` 为 `npm.cmd` 绝对路径（例如 `C:\Program Files\nodejs\npm.cmd`）后重试。

5. PowerShell 提示 `pyinstaller` 不是命令
- 使用模块方式执行：`python -m PyInstaller --onefile --name LabFlowLauncher labflow_launcher.py`。

6. EXE 提示找不到项目根目录
- 把 EXE 放到项目根目录，或设置 `LABFLOW_PROJECT_ROOT`。

7. 内网用户无法访问
- 检查服务器 IP、端口、防火墙规则。

8. 需要调试前端热更新
- 改用 `--frontend-mode dev`（只建议开发时）。
