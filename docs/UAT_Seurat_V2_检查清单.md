# UAT 检查清单：Seurat V2（Phase 4）

> 目标：使用至少 2 个不同 Seurat 数据集，验证 V2 流程（inspect → prepare-training → train 默认 prepared dataset）可用。

## 1. 测试前信息

- 测试日期：`____-__-__`
- 测试环境（本机/WSL/NAS+SSH）：`__________`
- 后端地址（默认）：`http://localhost:8000`
- 数据集 A（dataset_id）：`__________`
- 数据集 B（dataset_id）：`__________`
- `group_column`：`__________`
- `cluster_column`：`__________`
- `selected_clusters`：`__________`

## 2. 自动化 UAT 脚本（推荐先跑）

脚本路径：`scripts/uat_phase4_seurat_v2.py`

示例（完整链路，含训练）：

```bash
python scripts/uat_phase4_seurat_v2.py \
  --base-url http://localhost:8000 \
  --dataset-id <DATASET_ID_A> \
  --dataset-id <DATASET_ID_B> \
  --group-column sample \
  --cluster-column celltype \
  --selected-clusters T,B,NK \
  --seed 42
```

示例（仅检查 inspect + prepare，不跑训练）：

```bash
python scripts/uat_phase4_seurat_v2.py \
  --base-url http://localhost:8000 \
  --dataset-id <DATASET_ID_A> \
  --dataset-id <DATASET_ID_B> \
  --group-column sample \
  --cluster-column celltype \
  --selected-clusters T,B,NK \
  --seed 42 \
  --skip-train
```

脚本输出报告：`tmp_cache/uat_phase4_seurat_v2_report.json`

## 3. 手工检查清单（UI/API）

### 3.1 数据集 A

- [ ] `POST /api/seurat/inspect` 成功返回 metadata 与 UMAP（或返回无 UMAP 警告）
- [ ] `POST /api/seurat/prepare-training` 成功返回 `prepared_dataset_id`
- [ ] 预处理结果满足：`n_cells <= 500`、`n_genes <= 500`
- [ ] `POST /api/jobs/train` 提交成功
- [ ] 训练 job 中可见 `source_dataset_id`、`train_dataset_id`、`prepared_dataset_id`
- [ ] 若提供 `prepared_dataset_id`，训练实际使用对应 prepared 数据集

### 3.2 数据集 B

- [ ] `POST /api/seurat/inspect` 成功返回 metadata 与 UMAP（或返回无 UMAP 警告）
- [ ] `POST /api/seurat/prepare-training` 成功返回 `prepared_dataset_id`
- [ ] 预处理结果满足：`n_cells <= 500`、`n_genes <= 500`
- [ ] `POST /api/jobs/train` 提交成功
- [ ] 训练 job 中可见 `source_dataset_id`、`train_dataset_id`、`prepared_dataset_id`
- [ ] 若提供 `prepared_dataset_id`，训练实际使用对应 prepared 数据集

## 4. 验收标准（Phase 4）

- [ ] 至少 2 个数据集 UAT 通过
- [ ] `inspect`、`prepare-training`、`train` 三段链路可用
- [ ] 训练默认 prepared dataset 的行为符合预期
- [ ] 预处理摘要在前端可见且可追溯
- [ ] 报告与日志已归档（脚本 JSON + 训练 job log）

## 5. 结果记录

- 总结：`通过 / 不通过`
- 失败项与原因：`____________________`
- 修复任务链接：`____________________`
