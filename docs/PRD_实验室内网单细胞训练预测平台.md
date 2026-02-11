# Squidiff LabFlow（MVP）: 技术规范与 AI 指令

## 1. 需求审计总结（含缺失信息）

### 1.1 愿景（简化版）
- 这是一个实验室内网小工具，不是企业平台。
- 目标只有一件事：让几位实验室成员上传数据后，自动完成训练/预测并看到可视化结果。
- 不做多租户隔离、不做复杂权限体系、不做算法升级。

### 1.2 现有项目用途
- 当前 `Squidiff` 已具备脚本能力：读取 `h5ad`，训练模型（可选 `SMILES + dose`），并执行预测采样。

### 1.3 已训练模型后可以做什么
- 固定模型版本，给新样本做批量预测。
- 输出可视化报告（UMAP、热图、分组统计）供实验讨论。
- 在实验室内网共享使用，减少每个人手动跑脚本。

### 1.4 本期范围（MVP）
- 上传 Seurat 或 `h5ad`，可选上传 SMILES CSV。
- 自动做格式校验和转换（Seurat -> `h5ad`）。
- 一键训练或一键预测。
- 展示结果图并支持下载。

### 1.5 非目标
- 不做模型架构优化。
- 不做公网部署。
- 不做多团队数据隔离。
- 不做高并发和分布式训练。

### 1.6 缺失信息（需你后续确认）
- 服务器系统（Linux/Windows）和 GPU 规格。
- 单次数据规模（细胞数、基因数）上限。
- 是否需要最简单登录（共享账号）还是完全无登录。

## 2. 架构决策记录（含备选与权衡）

### 2.1 备选方案
- 方案 A：FastAPI + React + 文件存储（JSON）+ 本地任务队列（轻量）
- 优点：实现快、部署简单、足够支撑几个人使用。
- 缺点：扩展性一般，不适合大规模并发。
- 方案 B：FastAPI + Celery + Redis + 文件存储（JSON）+ React（标准化）
- 优点：更规范、更可扩展。
- 缺点：对当前“内网几个人”场景偏重。

### 2.2 选型结论（ADR-001）
- 选择方案 A（轻量优先）。
- 原则：先跑通稳定流程，再考虑重构成标准化多服务架构。

### 2.3 前端风格结论（ADR-002）
- 风格：科研报告感，清爽、克制、可读性优先。
- 配色：浅底 + 深色文字 + 单一强调色（蓝绿）。
- 页面：上传页、训练页、预测页、结果页四页即可，不做复杂门户。

## 3. 系统设计（目录结构 / 数据模型 / 关键流程）

### 3.1 目录结构（MVP）
```text
labflow/
  backend/
    app/
      api/
        datasets.py
        jobs.py
        results.py
      services/
        seurat_converter.py
        data_validator.py
        squidiff_runner.py
        visualize.py
      storage/
        state_manager.py
      main.py
    tests/
  frontend/
    src/
      pages/
      components/
      services/api.ts
      styles/tokens.css
  infra/
    docker-compose.yml
  docs/
    PRD_实验室内网单细胞训练预测平台.md
```

### 3.2 核心状态模型（最小，无 SQL）
- `datasets.json`：上传数据记录（原文件路径、转换后路径、校验状态）。
- `jobs.json`：训练/预测任务（状态、日志、参数、耗时）。
- `models.json`：模型输出（checkpoint 路径、关键参数）。
- `results.json`：预测结果与可视化文件路径。

### 3.3 API（最小可用）
- `POST /api/datasets/upload`：上传数据。
- `POST /api/datasets/{id}/validate`：校验字段与维度。
- `POST /api/jobs/train`：提交训练。
- `POST /api/jobs/predict`：提交预测。
- `GET /api/jobs/{id}`：查询状态与日志。
- `GET /api/results/{id}`：获取图和结果文件。

### 3.4 关键流程
- 上传：用户上传 Seurat/h5ad（可选 SMILES）-> 后端保存并校验。
- 训练：选择数据 -> 提交训练参数 -> 调用 `train_squidiff.py` -> 产出模型。
- 预测：选择模型和测试数据 -> 调用 `sample_squidiff.py` -> 输出预测矩阵。
- 可视化：自动生成 UMAP/热图等静态图并在前端展示。

## 4. 详细实现要求（错误处理 / 测试 / 安全 / 性能）

### 4.1 错误处理
- 所有任务状态仅四种：`queued/running/success/failed`。
- 失败时必须返回：失败原因 + 日志路径 + 重试按钮。
- 数据校验失败要明确提示字段缺失（如 `Group`、`SMILES`、`dose`）。

### 4.2 测试要求（MVP）
- 后端：至少覆盖上传校验、任务提交、状态查询。
- 前端：至少覆盖上传流程和任务状态展示。
- 全链路：至少 1 组小样本跑通“上传->训练->预测->可视化”。

### 4.3 安全要求（简化）
- 内网使用，默认可关闭复杂鉴权。
- 如需登录，先采用共享账号（单账号）即可。
- 上传文件限制类型和大小，防止误上传无关文件。

### 4.4 性能与部署基线（简化）
- 同时只跑 1 个训练任务；其余任务排队。
- 支持 3-5 人日常使用。
- 使用 `docker compose up -d` 单机部署。
- 状态持久化使用本地卷挂载（`state/*.json`），重启不丢任务记录。

## 5. AI 执行指令（分阶段任务清单）

### Phase 0: 工程初始化
- 创建 `backend/frontend/infra` 目录。
- 后端初始化 FastAPI；前端初始化 React + Vite。
- 定义最小 API 契约文档。
- Checkfix
- 1. `ruff check backend`
- 2. `ruff format --check backend`
- 3. `npm install`（依赖变更时）
- 4. `npm run lint --prefix frontend`
- 5. `npm run build --prefix frontend`

### Phase 1: 上传与校验
- 实现上传接口、文件落盘、Seurat 转换适配层。
- 实现 `Group` / `SMILES` / `dose` 校验。
- Checkfix
- 1. `ruff check backend`
- 2. `ruff format --check backend`
- 3. `uv run pytest -q backend/tests`

### Phase 2: 训练与预测封装
- 封装 `train_squidiff.py` 和 `sample_squidiff.py` 服务调用。
- 实现任务状态跟踪（排队、运行、完成、失败）。
- Checkfix
- 1. `ruff check backend`
- 2. `ruff format --check backend`
- 3. `uv run pytest -q backend/tests`

### Phase 3: 结果可视化与前端
- 前端实现 4 页：上传、训练、预测、结果。
- 结果页展示 UMAP/热图并支持下载结果。
- Checkfix
- 1. `npm run lint --prefix frontend`
- 2. `npm run build --prefix frontend`

### Phase 4: 内网部署与验收
- 编写 `docker-compose.yml` 并完成单机启动。
- 用真实小样本跑一次完整流程验收。
- Checkfix
- 1. `ruff check backend`
- 2. `ruff format --check backend`
- 3. `uv run pytest -q backend/tests`
- 4. `npm run lint --prefix frontend`
- 5. `npm run build --prefix frontend`

### MVP 验收标准（DoD）
- 实验室用户可在网页上完成：上传 -> 训练/预测 -> 查看图表 -> 下载结果。
- 失败任务可看到日志并可重试。
- 内网单机可稳定运行，重启后历史记录不丢失。
