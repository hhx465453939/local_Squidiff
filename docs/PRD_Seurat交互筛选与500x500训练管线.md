# Squidiff LabFlow V2: Seurat 交互筛选与 500x500 训练管线 PRD

## 1. 背景与目标

### 1.1 背景
- 当前 MVP 已支持上传 -> 校验 -> 训练 -> 结果查看，但对实验室用户最难的一步仍是“从 Seurat 对象中手工筛选合格细胞并准备训练输入”。
- 用户常见痛点：
- 1. 不清楚 metadata 应该怎么写。
- 2. 分组字段不固定（不一定叫 `Group`）。
- 3. 细胞和基因数量常超出训练可控范围。

### 1.2 V2 目标
- 在前端直接上传“已完成注释与 UMAP 的 Seurat 对象”后：
- 1. 自动渲染 UMAP；
- 2. 让用户指定哪个 metadata 列是分组字段；
- 3. 让用户选择要用于训练的细胞类群（cluster）；
- 4. 自动将训练输入裁剪为最多 `500 cells x 500 genes`；
- 5. 无缝进入下游训练、预测、可视化流程。

### 1.3 非目标
- 不改 Squidiff 模型架构与损失函数。
- 不引入数据库（继续 JSON 文件状态）。
- 不做公网多租户权限体系。

## 2. 用户输入与前置条件

### 2.1 必须输入
- `Seurat` 对象文件：`.rds` 或 `.h5seurat`。
- 对象必须已完成预处理与注释：
- 1. 有可用降维（至少 UMAP embeddings）；
- 2. 有细胞类群注释列（如 `celltype`, `cluster`, `cellclusters2`）；
- 3. 有实验分组注释列（如 `sample`, `condition`, `treatment`）。

### 2.2 可选输入
- `SMILES` 文件（CSV）用于药物结构模式。

### 2.3 Metadata 规范（用户文档）
- WebUI 中用户可指定：
- 1. `group_column`：分组字段（例如 `sample`）。
- 2. `cluster_column`：细胞类群字段（例如 `celltype`）。
- 后端会把 `group_column` 映射为训练所需 `Group`。
- 示例（R）：
```r
# 假设你的Seurat对象叫 seu
# 你的分组列是 sample，细胞类群列是 celltype
seu$Group <- seu$sample
table(seu$Group)
table(seu$celltype)
```

## 3. 核心流程（产品视角）

### 3.1 WebUI 流程
- Step A 上传 Seurat 对象。
- Step B 系统读取 metadata 列名与 UMAP。
- Step C 用户选择：
- 1. 分组字段（group_column）；
- 2. 细胞类群字段（cluster_column）；
- 3. 用于训练的目标细胞群集合（selected_clusters）。
- Step D 系统预处理并给出“拟训练矩阵尺寸预览”。
- Step E 用户确认后启动训练。
- Step F 训练完成后进入预测与可视化。

### 3.2 后端预处理流程
- 1. Seurat -> h5ad 转换。
- 2. 从 metadata 中抽取用户指定字段并标准化：
- `adata.obs["Group"] = adata.obs[group_column]`
- `adata.obs["Cluster"] = adata.obs[cluster_column]`
- 3. 过滤仅保留 `selected_clusters` 中细胞。
- 4. 细胞下采样（见算法 A）。
- 5. 特征基因筛选（见算法 B）。
- 6. 导出最终训练输入：最大 `500 x 500`。

## 4. 算法规格（必须实现）

### 4.1 算法 A：细胞抽样（最多 500）
- 输入：筛选后的细胞集合 `C`。
- 规则：
- 1. 若 `|C| <= 500`，全量保留。
- 2. 若 `|C| > 500`，执行“分层随机抽样（Stratified Sampling）”。
- 分层依据（优先顺序）：
- 1. `Group`（保证各分组比例）；
- 2. 若同组内仍过大，可按 `Cluster` 子分层。
- 抽样随机种子固定为可配置参数（默认 `42`），保证可复现。

### 4.2 算法 B：差异基因筛选（最多 500）
- 输入：抽样后细胞表达矩阵与 `Group`。
- 规则：
- 1. 若基因数 `<= 500`，全量保留。
- 2. 若基因数 `> 500`，按分组做差异分析，取显著性最高基因。
- 统计方法：
- 默认使用 Wilcoxon（与 Seurat 常规流程一致）。
- 产出：
- `top_genes = min(500, n_significant_genes)`。
- 最终训练矩阵上限：`500 cells x 500 genes`。

### 4.3 失败回退策略
- 若差异分析失败或显著基因不足：
- 回退到“高变基因（HVG）Top 500”。
- 若 UMAP 缺失：
- 不阻断训练流程，但前端提示“无 UMAP，跳过交互图”。

## 5. 系统设计变更

### 5.1 新增后端服务
- `backend/app/services/seurat_inspector.py`
- 读取 Seurat metadata 列、UMAP 坐标。
- `backend/app/services/dataset_preprocessor.py`
- 执行 cluster 过滤、分层抽样、差异基因筛选、导出训练 h5ad。

### 5.2 新增 API
- `POST /api/seurat/inspect`
- 输入：seurat 文件 ID。
- 输出：metadata 字段列表、UMAP 点数据预览、细胞总数。
- `POST /api/seurat/prepare-training`
- 输入：`dataset_id`, `group_column`, `cluster_column`, `selected_clusters`, `seed`。
- 输出：`prepared_dataset_id`, `n_cells`, `n_genes`, `sampling_report`, `gene_report`。
- `GET /api/seurat/prepare-training/{job_id}`
- 查看预处理任务状态与日志。

### 5.3 前端页面升级
- 上传页增加“Seurat 解析”状态。
- 新增“细胞筛选页”：
- 1. UMAP 交互图（框选/按类群勾选）；
- 2. 分组字段、类群字段选择器；
- 3. 预处理结果预览（最终矩阵尺寸）。
- 训练页默认使用 `prepared_dataset_id`。

## 6. 用户操作文档（给实验室成员）

### 6.1 你需要提前准备什么
- 一个预处理完成的 Seurat 对象：
- 1. 已有 UMAP；
- 2. 有细胞类群注释；
- 3. 有实验分组注释。

### 6.2 WebUI 实操步骤
- 1. 上传 `.rds`/`.h5seurat`。
- 2. 在“分组字段”里选择你的分组列（例如 `sample`）。
- 3. 在“细胞类群字段”里选择你的类群列（例如 `celltype`）。
- 4. 在 UMAP 或类群列表勾选要训练的细胞群。
- 5. 点击“准备训练数据”。
- 6. 确认显示 `n_cells <= 500` 且 `n_genes <= 500`。
- 7. 启动训练并等待结果。

### 6.3 常见错误
- “分组字段为空”：检查 metadata 列名拼写。
- “选定细胞太少”：建议每组至少 30 个细胞。
- “无 UMAP”：对象预处理阶段未保存降维结果。

## 7. 验收标准（DoD）

- 用户不需要写 R 代码也能完成细胞群选择并进入训练。
- 训练输入矩阵严格满足 `n_cells <= 500`、`n_genes <= 500`。
- 抽样与基因筛选报告可追踪（保存到 JSON + 日志）。
- 同一随机种子可复现实验结果（同输入得到同抽样集合）。
- 训练、预测、可视化全链路可跑通。

## 8. 实施计划（开发任务）

### Phase 1: Seurat 检查与 UI 选择
- 实现 `inspect` API。
- 前端增加 metadata 字段选择 + UMAP 展示。
- Checkfix:
- `ruff check backend`
- `ruff format --check backend`
- `npm run lint --prefix frontend`
- `npm run build --prefix frontend`

### Phase 2: 500x500 预处理管线
- 实现 `prepare-training` API。
- 落地分层抽样与差异基因筛选。
- 输出预处理报告（采样比例、入选基因清单）。
- Checkfix:
- `ruff check backend`
- `ruff format --check backend`
- `python -m compileall backend/app`

### Phase 3: 训练流程接入
- 训练任务默认接 `prepared_dataset_id`。
- 页面展示预处理摘要 + 训练来源可追溯信息。
- Checkfix:
- `ruff check backend`
- `ruff format --check backend`
- `npm run lint --prefix frontend`
- `npm run build --prefix frontend`

### Phase 4: 文档与实验室交付
- 补充 `docs/seurat转换指南.md` 的 V2 章节。
- 输出“实验室 10 分钟上手”文档。
- 小样本 UAT：至少 2 个不同 Seurat 数据集验证。
