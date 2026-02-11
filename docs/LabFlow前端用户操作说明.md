# LabFlow 前端用户操作说明

> 面向不熟悉命令行的用户，按页面步骤说明每个选项和参数如何填写。  
> 环境与启动方式见 `README.md` 与 `部署文档.md`。

---

## 流程概览

页顶流程：**上传 → 校验 → Seurat 检查 → 500x500 预处理 → 训练任务 → 结果页**。  
按顺序完成每一步即可；某步失败时，页面会显示错误信息，按提示修改后重试。

---

## 1) 上传数据

| 选项 | 说明 | 如何填写 |
|------|------|----------|
| **数据名称** | 给本数据集起的名字，便于区分 | 可选。例如：`colon_TC`、`筋膜_TC`。 |
| **主数据文件** | 单细胞数据文件 | 必选。支持 **.h5ad**、**.rds**、**.h5seurat**。点击「选择文件」选本地文件。 |
| **SMILES CSV（可选）** | 药物结构模式下的 SMILES + 剂量表 | 仅当勾选「使用药物结构模式」时需要。 |
| **使用药物结构模式** | 是否用 SMILES + dose 做药物结构训练 | 一般**不勾选**（基础转录组预测即可）。需要时再勾选并上传 SMILES CSV。 |

- 上传成功后，页面会显示 `dataset_id` 和 `status: uploaded`。
- 若主数据是 **.rds 或 .h5seurat**，下一步「校验」会触发 R 转 h5ad，需正确配置 R（见下一节）。

---

## 2) 校验

校验会检查数据格式，并在需要时把 .rds/.h5seurat 转成 .h5ad。**若 R 装在 Conda 里（常见于 Windows），必须按下面配置，否则会报「Rscript was not found」。**

### 选项说明

| 选项 | 说明 | 如何填写 |
|------|------|----------|
| **R 执行方式** | 后端如何调用 R | **Windows 且 R 在 Conda 中**：选 **cmd_conda**。<br>**Linux/mac 且 R 在系统 PATH 里**：可选 **direct**。 |
| **Rscript 命令** | direct 模式下要执行的命令 | direct 时多为 `Rscript`（或本机 Rscript 的完整路径）。cmd_conda 时可不改。 |
| **conda.bat 路径** | cmd_conda 时 Conda 的激活脚本 | **必须填完整路径**。例如：`F:\software\Miniconda3\condabin\conda.bat`。<br>可用 `where conda`（CMD）或 `Get-Command conda`（PowerShell）查本机路径；常见在 `...\Miniconda3\condabin\conda.bat` 或 `...\Anaconda3\condabin\conda.bat`。 |
| **Conda R 环境名** | cmd_conda 时要激活的 R 环境 | 填你本机 R 环境名，例如 **r-4.3**、**r-seurat**。可用 `conda info --envs` 查看。 |

### 推荐配置示例

- **Windows + Conda 安装的 R（如 r-4.3）**  
  - R 执行方式：**cmd_conda**  
  - conda.bat 路径：`F:\software\Miniconda3\condabin\conda.bat`（按本机实际路径改）  
  - Conda R 环境名：**r-4.3**

- **Linux/mac，R 已在系统 PATH**  
  - R 执行方式：**direct**  
  - Rscript 命令：**Rscript**

校验通过后会显示「校验结果: 通过」及 cells/genes 数量；若有错误，请根据红色提示修改 R 配置或数据后重试。

---

## 3) Seurat 解析（Seurat 检查）

- **作用**：读取数据的 metadata 字段和 UMAP 预览，供下一步「500x500 预处理」选列用。
- **前提**：已上传并**校验通过**（rds/h5seurat 会先转为 h5ad）。
- **操作**：点击「解析 Seurat」即可，无需额外参数。
- **结果**：展示 n_cells、n_genes、metadata 字段列表、UMAP 预览等；无 UMAP 不阻断后续步骤。

---

## 4) 500x500 预处理

在「准备训练数据（500x500）」前，需指定分组与类群列，以及要保留的类群。

| 选项 | 说明 | 如何填写 |
|------|------|----------|
| **分组字段（group_column）** | 用于分组的 metadata 列名 | 填列名，如 **Group**、**sample**、**condition**。可参考上一步解析出的「metadata 字段」。 |
| **类群字段（cluster_column）** | 细胞类型/聚类所在列 | 填列名，如 **Cluster**、**celltype**、**cell_type**。 |
| **随机种子（seed）** | 抽样与筛选的随机种子 | 一般用默认 **42** 即可，便于复现。 |
| **待筛选 clusters** | 要参与训练的类群标签，逗号分隔 | 填上一步解析出的该类群列中的**部分取值**，如 **T,B,NK**。留空或全选则通常表示全部类群参与（以实际后端逻辑为准）。 |

- 填写后点击「准备训练数据（500x500）」。
- 成功后会显示 **prepared_dataset_id**、n_cells、n_genes 等；后续训练会默认使用该 prepared 数据。

---

## 5) 提交训练任务

| 选项 | 说明 | 如何填写 |
|------|------|----------|
| **gene_size** | 输入基因数（与数据维度一致） | 校验通过后会自动带出建议值；也可手动改为与 prepared 数据一致。 |
| **output_dim** | 输出维度 | 一般与 **gene_size** 相同。 |
| **batch_size** | 训练批大小 | 默认 64；显存不足可改小（如 32）。 |
| **lr** | 学习率 | 默认 1e-4，一般可不改。 |

- 若上一步已完成预处理，页面会提示「当前训练默认使用 prepared_dataset_id: xxx」。
- 点击「开始训练」后，任务进入队列，可在「任务轮询状态」中查看进度。

---

## 6) 任务轮询状态

- 显示当前训练任务的 **job_id**、**status**（如 queued / running / success / failed）、**source_dataset_id**、**train_dataset_id**、**prepared_dataset_id**、开始/结束时间。
- 若 **status** 为 **failed**，会显示 **error_msg**，可根据提示排查（如 R 配置、数据维度、显存等）。
- 训练时间较长时请耐心等待；成功后进入「结果页」。

---

## 7) 结果页

- 训练 **status** 为 **success** 后，本区会显示：
  - **model_id**、**checkpoint** 路径、**gene_size**、**output_dim**；
  - 预测结果与图表（若有）；
  - 可展开「查看日志」查看训练日志。
- 图表与资产可通过页面提供的链接或后端 API 下载（见 API 文档）。

---

## 常见问题速查

| 现象 | 处理 |
|------|------|
| 校验报「Rscript was not found」 | 将「R 执行方式」改为 **cmd_conda**，并填写 **conda.bat 完整路径**和 **Conda R 环境名**（如 r-4.3）。 |
| 校验报「conda.bat 不是内部或外部命令」 | 后端未用最新代码或 conda.bat 路径错误。请填 conda.bat 的**完整路径**（如 F:\...\condabin\conda.bat），必要时重启后端。 |
| 预处理/训练报错 | 查看页面红色错误信息或「任务轮询状态」中的 error_msg；检查分组/类群列名、cluster 取值、gene_size/output_dim 是否与数据一致。 |
| 训练很慢或卡住 | 训练需 GPU；可查看后端日志或脚本输出。若为超时，见 README/全流程脚本中关于「训练轮询与 GPU 检测」的说明。 |

---

## 相关文档

- 环境安装与启动：`README.md` 第 5 节、`docs/部署文档.md`
- R/Seurat 转换细节：`docs/seurat转换指南.md`
- API 接口：`docs/api/seurat.md`、`backend/app/api/` 下各模块
