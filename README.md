# Squidiff LabFlow（本仓库说明）

> 这是一个“研究模型 + 内网 Web 工作流”的混合仓库。  
> 你既可以直接用 `train_squidiff.py` / `sample_squidiff.py` 做模型训练推理，也可以用前后端 Web 流程完成 Seurat 数据上传、500x500 预处理、训练与结果查看。

---

## 1. 项目功能总览

### 1.1 研究模型能力（根目录脚本）
- 基于扩散模型的单细胞转录组预测。
- 支持基础模式与药物结构模式（`SMILES + dose`）。
- 入口脚本：
  - `train_squidiff.py`
  - `sample_squidiff.py`

### 1.2 LabFlow Web 能力（`backend/` + `frontend/`）
- 数据上传与格式校验（支持 `.h5ad/.rds/.h5seurat`）。
- Seurat 检查（metadata 字段 + UMAP 预览）。
- 训练前预处理（Phase 2）：
  - cluster 过滤
  - 最多 500 cells 分层抽样
  - 最多 500 genes 筛选（Wilcoxon + fallback）
- 训练任务提交与轮询（Phase 3）：
  - 默认优先使用 `prepared_dataset_id`
  - 训练来源可追溯（`source_dataset_id`、`train_dataset_id`、`prepared_dataset_id`）
- 结果资产查看（模型信息、预测图像、日志）。

---

## 2. 当前架构（前后端 + 任务执行）

```text
Frontend (React/Vite)
   |
   | HTTP / JSON
   v
FastAPI backend
   ├─ /api/datasets   (上传/校验/转换)
   ├─ /api/seurat     (inspect + prepare-training)
   ├─ /api/jobs       (train/predict + poll + log)
   └─ /api/results    (模型/结果/资产)
   |
   v
JsonStateStore (backend/state/*.json)
   |
   v
JobQueue worker
   |
   v
SquidiffRunner -> train_squidiff.py / sample_squidiff.py
```

### 2.1 后端核心目录
- `backend/app/api/`：REST API 路由。
- `backend/app/services/`：业务服务层（转换、检查、预处理、任务执行）。
- `backend/app/storage/state_manager.py`：文件型状态存储。
- `backend/state/`：状态 JSON（`datasets/jobs/seurat_prepare_jobs/models/results`）。
- `backend/uploads/`：上传与预处理输出。
- `backend/artifacts/`：训练/预测任务产物与日志。

### 2.2 前端核心目录
- `frontend/src/App.tsx`：单页流程 UI（上传 -> 校验 -> inspect -> prepare -> train -> 结果）。
- `frontend/src/services/api.ts`：API 类型与请求封装。
- `frontend/src/styles/tokens.css`：样式 token 与页面样式。

---

## 3. API 概览

### 3.1 健康检查
- `GET /api/health`

### 3.2 数据集
- `GET /api/datasets`
- `POST /api/datasets/upload`
- `POST /api/datasets/{dataset_id}/validate`

### 3.3 Seurat（V2）
- `POST /api/seurat/inspect`
- `POST /api/seurat/prepare-training`
- `GET /api/seurat/prepare-training/{job_id}`

### 3.4 任务
- `GET /api/jobs`
- `GET /api/jobs/{job_id}`
- `GET /api/jobs/{job_id}/log`
- `POST /api/jobs/train`
- `POST /api/jobs/predict`

### 3.5 结果
- `GET /api/results`
- `GET /api/results/{result_id}`
- `GET /api/results/job/{job_id}`
- `GET /api/results/models/list`
- `GET /api/results/models/{model_id}`
- `GET /api/results/{result_id}/assets/{asset_name}`

详细接口请看：`docs/api/seurat.md`（Seurat 部分），其余接口可参考 `backend/app/api/*.py`。

---

## 4. 本地开发启动（推荐）

## 4.1 环境要求
- Python 3.11（推荐）
- Node.js 20+
- R + SeuratDisk（仅 `.rds/.h5seurat -> h5ad` 转换需要）

## 4.2 后端启动

```bash
pip install -r requirements.txt -r backend/requirements.txt
uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000
```

## 4.3 前端启动

```bash
cd frontend
npm install
npm run dev
```

默认访问：
- Frontend: `http://localhost:5173`
- Backend: `http://localhost:8000`

---

## 5. Docker 部署（内网快速拉起）

配置文件：`infra/docker-compose.yml`  
示例环境变量：`infra/.env.example`

```bash
cd infra
cp .env.example .env
docker compose up --build
```

### 5.1 关键环境变量

| 变量 | 说明 | 默认值 |
|---|---|---|
| `LABFLOW_DRY_RUN` | 后端是否走轻量 dry-run | `true` |
| `LABFLOW_R_EXEC_MODE` | R 执行模式（`direct`/`cmd_conda`） | `direct` |
| `LABFLOW_RSCRIPT_BIN` | Rscript 命令 | `Rscript` |
| `LABFLOW_R_CONDA_ENV` | `cmd_conda` 下的 Conda 环境名 | 空 |
| `LABFLOW_R_CONDA_BAT` | `conda.bat` 路径 | `conda.bat` |
| `VITE_API_BASE` | 前端请求后端地址 | `http://localhost:8000` |

---

## 6. 开发/运维常用命令

### 6.1 后端静态检查
```bash
ruff check backend/app backend/tests
ruff format --check backend/app backend/tests
```

### 6.2 前端检查与构建
```bash
cd frontend
npm run lint
npm run build
```

### 6.3 训练输入形状检查
```bash
python scripts/check_shape.py --data_path path/to/data.h5ad
```

### 6.4 Phase 4 UAT（至少 2 数据集）
```bash
python scripts/uat_phase4_seurat_v2.py \
  --base-url http://localhost:8000 \
  --dataset-id <A> \
  --dataset-id <B> \
  --group-column sample \
  --cluster-column celltype \
  --selected-clusters T,B,NK \
  --seed 42
```

---

## 7. 目录结构（简版）

```text
.
├─ backend/
│  ├─ app/
│  │  ├─ api/
│  │  ├─ services/
│  │  └─ storage/
│  ├─ state/
│  ├─ uploads/
│  ├─ artifacts/
│  └─ scripts/seurat_to_h5ad.R
├─ frontend/
│  └─ src/
├─ infra/
├─ docs/
├─ scripts/
├─ train_squidiff.py
└─ sample_squidiff.py
```

---

## 8. 文档导航

- 部署与环境：`docs/部署文档.md`
- Seurat 转换（含 V2 补充）：`docs/seurat转换指南.md`
- Seurat API：`docs/api/seurat.md`
- 10 分钟上手：`docs/实验室10分钟上手.md`
- UAT 清单：`docs/UAT_Seurat_V2_检查清单.md`
- 设计/需求：`docs/PRD_Seurat交互筛选与500x500训练管线.md`

---

## 9. 当前状态与注意事项

- V2 Phase 1~4 代码和交付文档已落地（inspect / prepare-training / train 默认 prepared dataset / UAT 资产）。
- 状态存储当前是 JSON 文件方案（MVP 取舍），不等同于数据库事务一致性。
- Seurat 转换依赖本机/容器内 R 运行时与 SeuratDisk 可用。
- 若你在本地开发，建议优先 `LABFLOW_DRY_RUN=true` 验证链路，再切真实训练。

---

## 10. 许可证与引用

- License: MIT（见 `LICENSE`）
- 论文引用见仓库根目录历史信息（`README` 旧版与论文条目）。
