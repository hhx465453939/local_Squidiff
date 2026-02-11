# Squidiff LabFlow（本仓库说明） / Squidiff LabFlow (Repository Guide)

<p align="center">
  <img src="docs/assets/labflow-logo.svg" alt="Squidiff LabFlow logo" width="760" />
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/siyuh/Squidiff/main/squidiff_logo.png" alt="Original Squidiff logo" width="190" />
</p>

> 致谢：本项目基于原始 Squidiff 工作持续扩展与工程化，感谢原作者及贡献者。  
> Acknowledgement: This project extends and engineers the original Squidiff work. Thanks to the original authors and contributors.
>
> 原项目地址：<https://github.com/siyuh/Squidiff>  
> Original project: <https://github.com/siyuh/Squidiff>

---

## 1. 项目定位 / Project Scope

本仓库是“研究脚本 + LabFlow Web 工作流”的混合仓库。  
This repository combines research scripts and the LabFlow web workflow.

你可以直接使用 `train_squidiff.py` / `sample_squidiff.py` 做训练和推理。  
You can directly use `train_squidiff.py` / `sample_squidiff.py` for training and inference.

你也可以通过前后端 Web 流程完成上传、校验、预处理、训练和结果查看。  
You can also complete upload, validation, preprocessing, training, and result review through the web flow.

---

## 2. 快速开始 / Quick Start

推荐顺序：上传 -> 校验 -> Seurat 检查 -> 500x500 预处理 -> 训练 -> 结果页。  
Recommended flow: upload -> validate -> Seurat inspect -> 500x500 prepare -> train -> results.

前端默认地址：`http://localhost:5173`。  
Frontend default URL: `http://localhost:5173`.

后端默认地址：`http://localhost:8000`。  
Backend default URL: `http://localhost:8000`.

前端用户操作说明见 `docs/LabFlow前端用户操作说明.md`。  
See `docs/LabFlow前端用户操作说明.md` for step-by-step frontend usage.

---

## 3. 核心能力 / Core Capabilities

研究脚本能力：扩散模型驱动的单细胞转录组预测。  
Research capability: diffusion-model-driven single-cell transcriptome prediction.

Web 能力：上传/校验、Seurat 检查、500x500 预处理、训练任务、结果资产查看。  
Web capability: upload/validation, Seurat inspect, 500x500 preprocessing, training jobs, and result assets.

任务调度支持按用户切换模式：`Serial (1)` 或 `Parallel (3)`。  
Task scheduling supports per-user mode switching: `Serial (1)` or `Parallel (3)`.

---

## 4. API 概览 / API Overview

健康检查：`GET /api/health`。  
Health check: `GET /api/health`.

数据集接口：`/api/datasets`。  
Dataset APIs: `/api/datasets`.

Seurat 接口：`/api/seurat`。  
Seurat APIs: `/api/seurat`.

任务接口：`/api/jobs`。  
Job APIs: `/api/jobs`.

结果接口：`/api/results`。  
Result APIs: `/api/results`.

用户调度偏好接口：`/api/user-prefs/scheduler`。  
User scheduler preference API: `/api/user-prefs/scheduler`.

---

## 5. 开发与部署（三选一） / Development and Deployment (Choose One)

环境可选：`uv`、`conda`、或 `venv + pip`。  
Environment options: `uv`, `conda`, or `venv + pip`.

先装 CUDA 版 PyTorch，再装 requirements，避免安装成 CPU 版本。  
Install CUDA PyTorch first, then install requirements, to avoid CPU-only fallback.

### 5.1 uv（推荐） / uv (Recommended)

```bash
# Windows
uv venv
.venv\Scripts\activate
uv pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
uv pip install -r requirements.txt -r backend/requirements.txt --extra-index-url https://download.pytorch.org/whl/cu118

# Linux / macOS
uv venv
source .venv/bin/activate
uv pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
uv pip install -r requirements.txt -r backend/requirements.txt --extra-index-url https://download.pytorch.org/whl/cu118
```

### 5.2 conda / conda

```bash
conda create -n labflow python=3.11
conda activate labflow
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install -r requirements.txt -r backend/requirements.txt --extra-index-url https://download.pytorch.org/whl/cu118
```

### 5.3 venv + pip / venv + pip

```bash
# Windows
python -m venv .venv
.venv\Scripts\activate
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install -r requirements.txt -r backend/requirements.txt --extra-index-url https://download.pytorch.org/whl/cu118

# Linux / macOS
python3 -m venv .venv
source .venv/bin/activate
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install -r requirements.txt -r backend/requirements.txt --extra-index-url https://download.pytorch.org/whl/cu118
```

### 5.4 启动后端 / Start Backend

```bash
python -m uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000
```

### 5.5 启动前端 / Start Frontend

```bash
cd frontend
npm install
npm run dev
```

---

## 6. 常用检查命令 / Common Check Commands

后端检查：  
Backend checks:

```bash
ruff check backend/app backend/tests
ruff format --check backend/app backend/tests
```

前端检查：  
Frontend checks:

```bash
cd frontend
npm run lint
npm run build
```

---

## 7. 文档导航 / Documentation Map

前端操作说明：`docs/LabFlow前端用户操作说明.md`。  
Frontend user guide: `docs/LabFlow前端用户操作说明.md`.

部署文档：`docs/部署文档.md`。  
Deployment guide: `docs/部署文档.md`.

Seurat 转换指南：`docs/seurat转换指南.md`。  
Seurat conversion guide: `docs/seurat转换指南.md`.

API 文档目录：`docs/api/`。  
API docs directory: `docs/api/`.

避坑指南：`docs/避坑指南.md`。  
Troubleshooting guide: `docs/避坑指南.md`.

---

## 8. 许可证 / License

本项目使用 MIT 许可证，详见 `LICENSE`。  
This project is released under the MIT License. See `LICENSE`.
