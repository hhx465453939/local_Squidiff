# CLAUDE.md — Squidiff / LabFlow 全栈开发指南

> 面向 AI 编码助手与开发者的项目上下文文档。  
> 本文重点覆盖：功能、架构、开发流程、部署、检查与交付。

---

## 1. 项目简介

`local_Squidiff` 是一个“研究模型 + 实验室 Web 工作流”融合仓库：

- **研究模型层**：扩散模型单细胞预测（`train_squidiff.py`, `sample_squidiff.py`）
- **应用层（LabFlow）**：前后端工作流，支持 Seurat 数据上传、校验、预处理、训练任务与结果展示

当前迭代重点是 LabFlow V2（Seurat 交互筛选 + 500x500 训练管线）。

---

## 2. 技术栈

### 2.1 后端
- Python + FastAPI
- Pydantic v2
- 文件型状态存储（JSON）
- 任务执行：后台线程队列 + 脚本调用

### 2.2 前端
- React 18 + TypeScript + Vite
- 单页流程 UI

### 2.3 数据/科学计算
- Scanpy / AnnData（读取和处理 `h5ad`）
- Seurat 转换脚本（R + SeuratDisk）

### 2.4 部署
- Docker + docker-compose（`infra/docker-compose.yml`）

---

## 3. 关键业务流程（V2）

1. 上传数据：`/api/datasets/upload`
2. 校验/转换：`/api/datasets/{dataset_id}/validate`
3. Seurat 检查：`/api/seurat/inspect`
4. 500x500 预处理：`/api/seurat/prepare-training`
5. 训练任务：`/api/jobs/train`（默认优先使用 prepared dataset）
6. 轮询任务：`/api/jobs/{job_id}`
7. 查看结果：`/api/results/*`

---

## 4. 后端架构

### 4.1 路由层（`backend/app/api/`）
- `datasets.py`：上传与校验
- `seurat.py`：inspect + prepare-training
- `jobs.py`：train/predict + log + 训练来源追踪
- `results.py`：结果、模型、资产文件访问

### 4.2 服务层（`backend/app/services/`）
- `data_validator.py`：数据校验
- `seurat_converter.py`：`rds/h5seurat -> h5ad`
- `seurat_inspector.py`：metadata/UMAP 检查
- `dataset_preprocessor.py`：cluster 过滤 + 分层抽样 + 基因筛选
- `job_queue.py`：任务队列
- `squidiff_runner.py`：训练/推理脚本执行

### 4.3 运行时与存储
- `runtime.py`：`store + runner + job_queue` 单例
- `storage/state_manager.py`：JSON CRUD + 时间戳

状态文件：
- `datasets.json`
- `jobs.json`
- `seurat_prepare_jobs.json`
- `models.json`
- `results.json`

---

## 5. 前端架构

### 5.1 页面组织
- 核心入口：`frontend/src/App.tsx`
- API 封装：`frontend/src/services/api.ts`
- 样式：`frontend/src/styles/tokens.css`

### 5.2 页面步骤（当前单页）
- 上传数据
- 校验
- Seurat 解析（字段 + UMAP）
- 500x500 预处理（字段/cluster/seed）
- 提交训练（默认带 `preparedDatasetId`）
- 任务状态与结果展示

---

## 6. 开发规范（必须遵守）

1. **API-first**：先后端 API，再前端消费。
2. **层级边界清晰**：业务逻辑不落在前端。
3. **最小改动原则**：不顺手改无关问题。
4. **文档同步**：API/流程变更后必须更新相关文档。
5. **记录调试**：更新 `.debug/labflow-mvp-debug.md`。

---

## 7. 本地开发运行

### 7.1 后端
```bash
pip install -r requirements.txt -r backend/requirements.txt
uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000
```

### 7.2 前端
```bash
cd frontend
npm install
npm run dev
```

---

## 8. Docker 部署

```bash
cd infra
cp .env.example .env
docker compose up --build
```

关键变量：
- `LABFLOW_DRY_RUN`
- `LABFLOW_R_EXEC_MODE`
- `LABFLOW_RSCRIPT_BIN`
- `LABFLOW_R_CONDA_ENV`
- `LABFLOW_R_CONDA_BAT`
- `VITE_API_BASE`

---

## 9. 检查与质量门禁

### 9.1 后端
```bash
ruff check backend/app backend/tests
ruff format --check backend/app backend/tests
```

### 9.2 前端
```bash
cd frontend
npm run lint
npm run build
```

### 9.3 UAT（Phase 4）
```bash
python scripts/uat_phase4_seurat_v2.py --help
```

---

## 10. 文档索引

- 总览：`README.md`
- 部署：`docs/部署文档.md`
- Seurat 转换：`docs/seurat转换指南.md`
- Seurat API：`docs/api/seurat.md`
- 10 分钟上手：`docs/实验室10分钟上手.md`
- UAT 清单：`docs/UAT_Seurat_V2_检查清单.md`
- 调试记录：`.debug/labflow-mvp-debug.md`

---

## 11. 当前已知限制

- JSON 文件存储适合 MVP，不适合高并发生产写入。
- Seurat 转换依赖 R 运行时可用性。
- 完整训练链路受本机模型依赖与硬件资源影响较大（可先用 `LABFLOW_DRY_RUN=true` 验证流程）。

---

## 12. 给 AI 代理的执行建议

- 先判定问题归属层，再动手改代码。
- 多层需求按顺序执行：后端 -> 文档 -> 前端 -> 集成。
- 修改 API 契约后，必须同步：
  - `docs/api/*`
  - `README.md`（若影响用户路径）
  - `.debug/labflow-mvp-debug.md`
