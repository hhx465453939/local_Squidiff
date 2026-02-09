# AGENTS.md - local_Squidiff 协作与开发规范

> 作用：给 AI Coding Agent 与协作者提供统一上下文。  
> 范围：本文件所在目录及其子目录。  
> 当前仓库路径：`E:\Development\local_Squidiff`

---

## 1. 项目定位

本仓库包含两条并行能力线：

1. **研究脚本线（根目录）**  
   - `train_squidiff.py` / `sample_squidiff.py`
   - 面向模型训练与推理

2. **LabFlow Web 线（前后端）**  
   - `backend/` + `frontend/` + `infra/`
   - 面向实验室内网工作流（上传、校验、Seurat 检查、500x500 预处理、训练任务、结果展示）

---

## 2. 架构与层级边界（强约束）

### 2.1 三层职责
- **前端层 (`frontend/`)**：页面渲染、用户交互、调用 API、状态展示。
- **后端层 (`backend/app/`)**：业务逻辑、数据处理、任务调度、API 暴露。
- **存储执行层 (`backend/state/`, `backend/uploads/`, `backend/artifacts/`)**：
  - JSON 状态存储
  - 上传数据、预处理输出
  - 训练与预测产物

### 2.2 禁止事项
- 不在前端实现后端业务逻辑。
- 不用前端绕过后端 bug。
- 不跨层做临时 workaround 代替根因修复。

---

## 3. API-First 开发流程（强约束）

后端每个功能按以下 5 步闭环执行：

1. Implement（实现）
2. Checkfix（lint/format/build）
3. Encapsulate（服务层封装）
4. Expose API（路由与契约）
5. Document API（同步文档）

跨层任务执行顺序固定：

1. 后端实现
2. API 文档更新
3. 前端接入
4. 集成验证

---

## 4. 运行与部署上下文

### 4.1 本地开发（默认）
- Backend: `uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000`
- Frontend: `cd frontend && npm run dev`

### 4.2 Docker 部署
- 配置文件：`infra/docker-compose.yml`
- 环境变量模板：`infra/.env.example`

### 4.3 R 转换执行模式
- `LABFLOW_R_EXEC_MODE=direct|cmd_conda`
- `cmd_conda` 需要正确设置：
  - `LABFLOW_R_CONDA_ENV`
  - `LABFLOW_R_CONDA_BAT`

---

## 5. 当前后端模块与 API（速查）

### 5.1 `backend/app/api/datasets.py`
- 上传与校验
- `POST /api/datasets/upload`
- `POST /api/datasets/{dataset_id}/validate`

### 5.2 `backend/app/api/seurat.py`
- V2 inspect + prepare-training
- `POST /api/seurat/inspect`
- `POST /api/seurat/prepare-training`
- `GET /api/seurat/prepare-training/{job_id}`

### 5.3 `backend/app/api/jobs.py`
- 训练/预测任务
- `POST /api/jobs/train`（默认优先使用 prepared dataset）
- `POST /api/jobs/predict`
- `GET /api/jobs/{job_id}`
- `GET /api/jobs/{job_id}/log`

### 5.4 `backend/app/api/results.py`
- 结果、模型、资产访问

---

## 6. 状态与数据文件（不可忽略）

- `backend/state/datasets.json`
- `backend/state/jobs.json`
- `backend/state/seurat_prepare_jobs.json`
- `backend/state/models.json`
- `backend/state/results.json`

上传/产物目录：
- `backend/uploads/`
- `backend/artifacts/`

---

## 7. 文档同步规则（强约束）

当你修改以下内容时，必须同步文档：

- API 变更 -> 更新 `docs/api/*.md`
- 数据流程或部署变更 -> 更新 `README.md` 与 `docs/部署文档.md`
- Seurat 流程变更 -> 更新 `docs/seurat转换指南.md`
- 交付与验收变更 -> 更新 `docs/UAT_Seurat_V2_检查清单.md`

---

## 8. Checkfix 规则

### 8.1 Python
```bash
ruff check backend/app backend/tests
ruff format --check backend/app backend/tests
```

### 8.2 Frontend
```bash
cd frontend
npm run lint
npm run build
```

### 8.3 UAT 脚本（Phase 4）
```bash
python scripts/uat_phase4_seurat_v2.py --help
```

---

## 9. Skills（可用能力）

本仓库内置 skills 目录：`.codex/skills/`

| Skill | 场景 | 路径 |
|---|---|---|
| `api-first-modular` | 前后端分层、API 优先、跨层分解 | `.codex/skills/api-first-modular/` |
| `code-debugger` | Bug 修复、增量开发、`.debug/` 记录 | `.codex/skills/code-debugger/` |
| `debug-ui` | 前端样式与交互调试 | `.codex/skills/debug-ui/` |
| `ai-spec` | 需求转技术规格 | `.codex/skills/ai-spec/` |
| `ralph` | PRD 驱动任务循环 | `.codex/skills/ralph/` |

---

## 10. Skill 触发规则

- 用户显式提到 skill 名（如 `$code-debugger`）-> 必须使用该 skill。
- 任务语义明显匹配 skill 描述 -> 应主动使用。
- 同时命中多个 skill -> 选择最小必要集合并说明执行顺序。

---

## 11. `.debug/` 记录规则

- 按模块维护调试历史，统一放在 `.debug/`。
- 每次非微小改动应更新对应 debug 文档，记录：
  - 问题与根因
  - 变更文件与函数
  - 检查命令与结果
  - 影响面与后续事项

当前主记录：`.debug/labflow-mvp-debug.md`
