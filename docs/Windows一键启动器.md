# Windows 一键启动器（LabFlow 实验室内网）

> 目标：给不会敲命令行的实验室同事一个「双击即可用」的入口，在本机启动 LabFlow，并允许内网其他人通过浏览器访问和注册账号。

---

## 1. 前置条件

在使用一键启动器之前，请先完成标准环境准备（只需要 **技术同学做一次**）：

- 操作系统：Windows 10 或以上  
- 已克隆本仓库到本机（例如 `E:\Development\local_Squidiff`）  
- 已按 `README` 中的「5. 开发与部署」完成以下步骤：
  - 准备 Python 环境（`uv` / `conda` / venv 三选一）
  - 安装后端依赖：
    - `pip install -r requirements.txt -r backend/requirements.txt`
  - 在 `frontend/` 目录下至少执行过一次：
    - `npm install`

> 启动器只是**帮你自动敲命令**，并不会替代依赖安装流程。如果依赖未安装好，启动器会直接报错退出。

---

## 2. 直接用 Python 启动（开发者/管理员）

适合技术同学在自己的机器上先验证：

1. 打开 PowerShell 或 CMD  
2. 进入项目根目录（根据你的路径调整）：

   ```powershell
   cd C:\Users\DamnCheater\.cursor\worktrees\local_Squidiff\cpx
   ```

3. 激活你为 LabFlow 准备的 Python 环境（任选其一）：

   - **uv 环境**：

     ```powershell
     .venv\Scripts\activate
     ```

   - **conda 环境**：

     ```powershell
     conda activate labflow
     ```

   - **本机 venv**：

     ```powershell
     .venv\Scripts\activate
     ```

4. 运行启动器脚本：

   ```powershell
   python labflow_launcher.py
   ```

5. 等待控制台中出现 Uvicorn 和 Vite 启动日志后，你会看到类似提示：

   - 本机访问：`http://localhost:5173`、`http://localhost:8000`
   - 内网访问：`http://<本机IP>:5173`、`http://<本机IP>:8000`

6. 要停止服务时，在该窗口按下 `Ctrl + C`，等待提示「LabFlow 服务已停止」后再关闭窗口。

---

## 3. 打包为 .exe（给非技术同事分发）

如果你希望给同事一个 **双击即可运行** 的 `LabFlowLauncher.exe`，可以使用 PyInstaller 将 `labflow_launcher.py` 打包成单文件可执行程序。

> 注意：打包后的 .exe **只负责帮你启动现有 Python/Node 环境**，它不会打包整个 LabFlow 项目代码和依赖。  
> 因此：同事的机器上仍然需要事先部署好 LabFlow（推荐在一台「服务器」上运行，然后通过内网访问）。

### 3.1 安装 PyInstaller

在你用于运行 LabFlow 的同一个 Python 环境中安装 PyInstaller：

```powershell
pip install pyinstaller
```

### 3.2 打包启动器

在项目根目录执行：

```powershell
pyinstaller --onefile --name LabFlowLauncher labflow_launcher.py
```

成功后会在 `dist/` 目录下生成：

- `dist/LabFlowLauncher.exe`

你可以把这个文件复制到桌面、U 盘或共享文件夹中发给同事使用。

### 3.3 使用 .exe 启动

假设你已经在某台 Windows 机器上按 README 完整部署好 LabFlow，并且该机器用于实验室内网共享使用：

1. 将 `LabFlowLauncher.exe` 放在 **项目根目录**（与 `labflow_launcher.py` 同级）或任意位置  
2. 双击运行 `LabFlowLauncher.exe`  
3. 等待窗口中出现访问地址提示：
   - 本机浏览器输入：`http://localhost:5173`
   - 内网其他机器输入：`http://<这台机器的 IP>:5173`
4. 要停止服务时，关闭窗口或在窗口中按 `Ctrl + C`。

> 建议：在用于长期服务的机器上，将 `LabFlowLauncher.exe` 固定在任务栏或创建桌面快捷方式，方便实验室成员一键启动。

---

## 4. 内网访问与注册使用方式

启动器使用的启动命令与 README 完全一致：

- 后端：`python -m uvicorn backend.app.main:app --host 0.0.0.0 --port 8000`
- 前端（Vite）：监听 `0.0.0.0:5173`（见 `frontend/vite.config.ts`）

因此，只要：

- 这台 Windows 机器处于实验室内网中；
- 防火墙允许 5173 和 8000 端口被局域网访问；

那么其他同事即可通过：  

- 前端：`http://<这台机器 IP>:5173`  
- 后端 API（一般由前端自动调用）：`http://<这台机器 IP>:8000`

进行访问。

### 4.1 注册与登录

LabFlow 已内置账号系统（SQLite 本地数据库 + 加密密码）：

- 首次使用时，在首页直接 **注册新账号**；
- 注册后使用该账号登录，点击「开始分析」即可进入完整流程；
- 登录后前端会自动携带 Bearer Token 调用后端 API。

如需查看账号系统相关 API 文档，请参考：

- `docs/api/auth.md`

---

## 5. 常见问题与排查

| 问题 | 可能原因 | 排查与解决 |
|------|----------|------------|
| 启动器提示找不到 `python` | Python 未加入 PATH，或使用了特殊路径 | 在当前 PowerShell/CMD 中能否直接运行 `python --version`；如不行，请在启动器前先激活虚拟环境，或设置环境变量 `LABFLOW_PYTHON` 为完整 Python 路径 |
| 启动器提示找不到 `npm` | Node.js 未安装或未加入 PATH | 安装 Node.js 20+，重新打开终端，确认 `npm -v` 正常输出 |
| 浏览器访问 `http://<IP>:5173` 打不开 | Windows 防火墙或安全软件拦截端口 | 临时关闭防火墙测试，或为 5173/8000 开放入站规则 |
| 前端页面打开但报「无法连接后端」 | 后端未成功启动，或端口/主机不对 | 查看启动器窗口中 Uvicorn 日志是否有报错；确认 8000 端口未被占用 |
| 登录/注册失败 | 数据库文件权限或路径异常 | 检查后端日志，确认 SQLite 数据库文件可读写；如需迁移数据库路径，请参考后端配置 |

---

## 6. 建议使用方式

- 开发调试：技术同学直接运行 `python labflow_launcher.py`，享受一键启动 + 自动打印本机/内网地址。  
- 实验室日常使用：在一台固定 Windows 机器上部署 LabFlow + `LabFlowLauncher.exe`，把该机器当作「小服务器」，同事只需记住一个内网地址即可使用并自行注册账号。

