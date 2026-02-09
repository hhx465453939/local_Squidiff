# Seurat 转换指南（R -> h5ad）

> 目标：将 Seurat 对象转换为 Squidiff 可直接使用的 `h5ad`，并完成最小验证。

## 1. 必要数据要求

转换后数据需满足：

- `adata.X` 存在，且维度为 `(n_cells, n_genes)`
- `adata.obs["Group"]` 存在
- 若使用药物结构：`adata.obs["SMILES"]` 与 `adata.obs["dose"]` 必须存在

训练参数约束：

- `--gene_size == adata.n_vars`
- `--output_dim == adata.n_vars`

## 2. 推荐转换路径（SeuratDisk）

### 2.1 R 端导出

```r
library(Seurat)
library(SeuratDisk)

# 确保 Group 存在
if (!"Group" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$Group <- seurat_obj$seurat_clusters
}

# 推荐流程：先 h5seurat，再转 h5ad
SaveH5Seurat(seurat_obj, filename = "sample.h5seurat", overwrite = TRUE)
Convert("sample.h5seurat", dest = "h5ad", overwrite = TRUE)
```

### 2.2 Python 端验证

```python
import scanpy as sc

adata = sc.read_h5ad("sample.h5ad")
print("cells:", adata.n_obs)
print("genes:", adata.n_vars)
print("obs columns:", adata.obs.columns.tolist())

assert "Group" in adata.obs.columns, "missing Group"
```

## 3. 药物结构场景补充

当训练使用 `--use_drug_structure True` 时：

1. 处理组数据中必须有 `SMILES` 和 `dose`。  
2. 需要提供 `--control_data_path`。  
3. 对照组与处理组的基因维度必须一致。  

SMILES 简单校验示例：

```python
from rdkit import Chem

def valid_smiles(s):
    return Chem.MolFromSmiles(s) is not None
```

## 4. 常见问题

### 4.1 `KeyError: 'Group'`

- 在 R 中补齐 `Group` 列后重新导出。

### 4.2 维度不匹配

- 用 `adata.n_vars` 回填训练命令里的 `--gene_size` 与 `--output_dim`。

### 4.3 `SMILES` 无效

- 用 RDKit 先批量过滤非法 SMILES，再训练。

## 5. 转换后检查清单

- [ ] `adata.n_obs > 0`
- [ ] `adata.n_vars > 0`
- [ ] `Group` 列存在
- [ ] 若使用药物结构：`SMILES` 与 `dose` 列存在
- [ ] `--gene_size` 与 `--output_dim` 已和 `adata.n_vars` 对齐
- [ ] 训练前可通过 `python scripts/check_shape.py --data_path ...` 复核

## 6. V2（Seurat 交互筛选 + 500x500 训练管线）补充

> 适用于 `docs/PRD_Seurat交互筛选与500x500训练管线.md` 的 Phase 1~3。

### 6.1 V2 数据要求（在原有要求基础上新增）

- 推荐在 Seurat 中提前包含：
  - 分组字段（例如 `sample`，用于 `group_column`）
  - 类群字段（例如 `celltype`，用于 `cluster_column`）
  - UMAP（用于前端交互预览）
- 若没有 UMAP：不阻断后续训练，但前端会提示“无 UMAP，跳过交互图”。

### 6.2 建议的 metadata 规范化（R 侧）

```r
# 建议统一列名，减少前端选择歧义
if (!"Group" %in% colnames(seurat_obj@meta.data) && "sample" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$Group <- seurat_obj$sample
}
if (!"Cluster" %in% colnames(seurat_obj@meta.data) && "celltype" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$Cluster <- seurat_obj$celltype
}
```

### 6.3 V2 后端流程（接口顺序）

1. 上传并校验（必要时自动转换为 `h5ad`）  
2. `POST /api/seurat/inspect`：读取 metadata + UMAP 预览  
3. `POST /api/seurat/prepare-training`：cluster 过滤 + 500 cells 抽样 + 500 genes 筛选  
4. `POST /api/jobs/train`：默认优先使用 `prepared_dataset_id`  

参考接口文档：`docs/api/seurat.md`。

### 6.4 500x500 约束确认

- 预处理结果需满足：
  - `n_cells <= 500`
  - `n_genes <= 500`
- 若使用前端流程，训练参数 `gene_size` 与 `output_dim` 会自动回填为预处理后的 `n_genes`。

### 6.5 API 快速自检（可选）

```bash
curl -X POST http://localhost:8000/api/seurat/inspect \
  -H "Content-Type: application/json" \
  -d '{"dataset_id":"<dataset-id>"}'
```

```bash
curl -X POST http://localhost:8000/api/seurat/prepare-training \
  -H "Content-Type: application/json" \
  -d '{
    "dataset_id":"<dataset-id>",
    "group_column":"sample",
    "cluster_column":"celltype",
    "selected_clusters":["T","B"],
    "seed":42
  }'
```
