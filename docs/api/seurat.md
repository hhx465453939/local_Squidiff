# Seurat API 文档

> 最后更新: 2026-02-10 | 版本: v1.0

## 概览
Seurat API 提供两类能力：
- 数据检查：读取 metadata 字段与 UMAP 预览；
- 训练预处理：按用户选择执行 cluster 过滤、500 细胞分层抽样、500 基因筛选，并输出可训练 `h5ad` 与报告。

## Base URL
`/api/seurat`

## 认证方式
无认证（当前 MVP 内网模式）。

## 端点列表
| 方法 | 路径 | 功能 | 认证 |
|------|------|------|------|
| POST | `/api/seurat/inspect` | 检查指定数据集的 metadata 与 UMAP 预览 | 否 |
| POST | `/api/seurat/prepare-training` | 执行 500x500 训练预处理并生成新数据集 | 否 |
| GET | `/api/seurat/prepare-training/{job_id}` | 查询预处理任务状态与结果 | 否 |

## 详细接口

### 检查 Seurat 数据集
- **路径**: `POST /api/seurat/inspect`
- **功能**: 从已准备好的 `h5ad` 读取基础信息、metadata 列名和 UMAP 预览点。
- **认证**: 否

#### 请求参数
| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| `dataset_id` | string | 是 | - | 数据集 ID（来自 `/api/datasets/upload`） |
| `umap_preview_limit` | integer | 否 | 1500 | UMAP 预览点上限，范围 `[1, 5000]` |

#### 响应格式
**成功 (200)**:
```json
{
  "dataset_id": "uuid",
  "inspect": {
    "n_cells": 1234,
    "n_genes": 5678,
    "metadata_columns": ["Group", "celltype"],
    "has_umap": true,
    "umap": {
      "key": "X_umap",
      "n_points": 1234,
      "preview_count": 500,
      "truncated": true,
      "preview": [
        { "cell_id": "AAAC...", "x": 1.23, "y": -0.45 }
      ]
    },
    "warnings": []
  }
}
```

**失败 (4xx)**:
```json
{ "detail": "Dataset not found" }
```

#### 错误说明
| HTTP状态码 | 含义 | 处理建议 |
|-----------|------|----------|
| 404 | 数据集不存在 | 检查 `dataset_id` 是否正确 |
| 400 | 缺少 `path_h5ad` 或文件不存在 | 先调用 `/api/datasets/{dataset_id}/validate` 触发转换 |
| 400 | 缺少运行依赖（如 `scanpy`） | 在后端环境安装依赖后重试 |

#### 调用示例
```bash
curl -X POST http://localhost:8000/api/seurat/inspect \
  -H "Content-Type: application/json" \
  -d '{"dataset_id":"<dataset-id>","umap_preview_limit":500}'
```

### 执行训练预处理
- **路径**: `POST /api/seurat/prepare-training`
- **功能**:
  1) 按 `selected_clusters` 过滤细胞；
  2) 按 `Group -> Cluster` 分层抽样到最多 500 cells；
  3) 优先用 Wilcoxon DEG 筛选最多 500 genes，失败回退 HVG/方差法；
  4) 产出预处理 `h5ad`，注册为新 `prepared_dataset_id`。
- **认证**: 否

#### 请求参数
| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| `dataset_id` | string | 是 | - | 原始数据集 ID（需已有 `path_h5ad`） |
| `group_column` | string | 是 | - | 分组字段列名 |
| `cluster_column` | string | 是 | - | 类群字段列名 |
| `selected_clusters` | string[] | 是 | - | 要保留的 cluster 值列表 |
| `seed` | integer | 否 | 42 | 随机种子，保证抽样可复现 |

#### 响应格式
**成功 (200)**:
```json
{
  "job_id": "uuid",
  "prepared_dataset_id": "uuid",
  "n_cells": 500,
  "n_genes": 500,
  "sampling_report": {
    "mode": "stratified_sampling",
    "seed": 42
  },
  "gene_report": {
    "method": "wilcoxon_deg",
    "fallback_used": false
  }
}
```

**失败 (4xx)**:
```json
{ "detail": "No cells left after cluster filtering" }
```

#### 错误说明
| HTTP状态码 | 含义 | 处理建议 |
|-----------|------|----------|
| 404 | 数据集不存在 | 检查 `dataset_id` |
| 400 | 缺少 `path_h5ad` | 先执行 `/api/datasets/{dataset_id}/validate` |
| 400 | 字段不存在/过滤后无细胞 | 检查 `group_column`、`cluster_column`、`selected_clusters` |
| 400 | 依赖缺失（如 `scanpy`） | 安装依赖后重试 |

#### 调用示例
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

### 查询预处理任务
- **路径**: `GET /api/seurat/prepare-training/{job_id}`
- **功能**: 返回任务状态、错误信息（若失败）和预处理结果（若成功）。
- **认证**: 否

#### 响应格式
**成功 (200)**:
```json
{
  "job": {
    "id": "uuid",
    "status": "success",
    "prepared_dataset_id": "uuid",
    "result": {
      "n_cells": 500,
      "n_genes": 500
    }
  }
}
```
