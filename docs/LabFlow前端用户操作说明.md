# LabFlow 前端用户操作说明

本文档面向实验室内网用户，按“从登录到训练结果”的顺序给出可直接照做的步骤。

---

## 0. 启动前置条件

请先确认后端与前端都已启动。

1. 启动后端（项目根目录）
```bash
python -m uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000
```

2. 启动前端（新终端）
```bash
cd frontend
npm run dev
```

3. 打开前端页面（通常是）
- `http://localhost:5173`

4. 健康检查
- 浏览器可访问：`http://localhost:8000/api/health`
- 期望返回：`{"status":"ok"}`

如果出现 `ERR_CONNECTION_REFUSED`，说明后端没启动或端口不是 8000。

---

## 1. 首页登录与入口

### 1.1 注册并登录

1. 打开首页后，在“账户入口”输入用户名和密码。
2. 首次使用点击“注册并登录”。
3. 后续使用点击“登录”。

账号规则：
- 用户名：3~32 位，只能包含字母、数字、下划线、点、短横线。
- 密码：至少 8 位。

### 1.2 登录后行为

- 登录成功后，页面会自动进入分析流程页（无需再点“开始分析”）。
- 首页与分析页都保留“用户说明书”入口。
- 左侧有 `Task Center` 侧边栏，可管理任务/模型/报告。

### 1.3 Task Center（侧边栏）

侧边栏提供三个区块：

1. `Jobs`
- 查看历史任务与正在运行任务。
- `Serial (1) / Parallel (3)`：切换“当前用户”的任务调度模式。
  - `Serial (1)`：同一时间最多 1 个任务运行。
  - `Parallel (3)`：同一时间最多 3 个任务运行。
- 新账号默认是 `Parallel (3)`，可随时切回 `Serial (1)`。
- `My scheduler: ...`：显示当前模式与当前用户运行中任务数量（例如 `parallel (2/3)`）。
- `Open`：恢复该任务的轮询、日志与状态显示。
- `Delete`：删除已结束任务（运行中任务需先停止）。
- `New task`：清空当前表单，开始新任务流程。
- `One-click clear`：一键清理卡住的 `queued/running` 测试任务（会连带删除相关产物）。
- `Refresh`：立即拉取最新任务列表。
- `Logout`：退出当前账号，返回首页，可切换其他账号登录。

2. `Models`
- 查看已产出的模型记录。
- `Job`：跳转并恢复关联任务视图。
- `Download`：下载模型 checkpoint。
- `Delete`：删除模型记录与文件。

3. `Reports`
- 查看结果报告记录。
- `Job`：跳转并恢复关联任务视图。
- `Download`：下载报告资产（优先 `summary.json`）。
- `Delete`：删除结果记录与报告文件。

---

## 2. 上传数据

在“1) 上传数据”中完成：

1. （可选）填写数据名称。
2. 上传主数据文件（支持 `.h5ad` / `.rds` / `.h5seurat`）。
3. （可选）上传 SMILES CSV。
4. 按需勾选“使用药物结构模式”。
5. 点击“上传”。

上传后会显示 `dataset_id` 与 `status`。

---

## 3. 数据校验

在“2) 校验”中完成运行参数设置后点击“开始校验”。

### 3.1 Windows + Conda（推荐）

建议配置：
- R 执行方式：`cmd_conda`
- conda.bat：完整路径（示例）`F:\software\Miniconda3\condabin\conda.bat`
- Conda R 环境名：示例 `r-4.3`

### 3.2 direct 模式

当系统 PATH 已可直接调用 `Rscript` 时，可选：
- R 执行方式：`direct`
- Rscript 命令：`Rscript`

---

## 4. Seurat 检查（V2）

在“3) Seurat 解析”点击“解析 Seurat”。

可查看：
- `n_cells` / `n_genes`
- metadata 字段
- UMAP 预览（若数据包含）

点击 metadata 字段后，前端会显示取值及每类细胞计数（类似 `table()`）。

---

## 5. 500x500 预处理

在“4) 500x500 预处理”填写：
- `group_column`
- `cluster_column`
- `seed`
- 待筛选 clusters（逗号分隔）

点击“准备训练数据（500x500）”后，成功会得到：
- `prepared_dataset_id`
- 预处理后细胞数与基因数

---

## 6. 提交训练任务

在“5) 提交训练任务”设置：
- `gene_size`
- `output_dim`
- `batch_size`
- `lr`

点击“开始训练”。

---

## 7. 任务轮询、日志和 GPU 看板

在“6) 任务轮询状态”：

1. 查看 `job_id`、状态、开始/结束时间。
2. 任务 `running` 时可点击“停止训练任务”。
3. 运行日志会实时刷新（无需先停止任务）。
4. 若机器有 NVIDIA GPU，会显示 `GPU Runtime`：
   - GPU 利用率
   - 显存占用
   - 显存利用率
   - 温度
   - 功耗

如果没有安装 `nvidia-smi` 或机器无 NVIDIA，GPU 面板会显示不可用原因，但训练流程不受影响。

---

## 8. 结果查看

训练成功后，在“7) 结果页”可查看：
- 模型信息（model_id、checkpoint、参数）
- 结果资产（图像等）
- 任务日志

---

## 9. 常见问题

### 9.1 页面提示 `ERR_CONNECTION_REFUSED`

原因：后端未启动或端口不对。

处理：
1. 重新启动后端。
2. 确认访问 `http://localhost:8000/api/health` 返回 `{"status":"ok"}`。

### 9.2 登录/注册失败

1. 检查用户名和密码是否满足规则。
2. 查看后端终端是否有报错。
3. 若报 401，通常是 token 失效，重新登录即可。

### 9.5 任务显示 running 但无法停止/删除

通常是历史测试任务在后端重启后变成“假 running”状态。

处理：
1. 打开左侧 `Task Center`。
2. 点击 `One-click clear`。
3. 等待清理完成后点击 `Refresh`。
4. 再重新提交新任务。

### 9.3 用户说明书打不开

1. 确认后端已启动。
2. 直接访问：`http://localhost:8000/api/auth/user-guide`
3. 若需原始 Markdown：`http://localhost:8000/api/auth/user-guide?raw=true`

### 9.4 日志长时间显示“等待日志”

1. 先确认任务状态是否已进入 `running`。
2. 检查后端目录：`backend/artifacts/jobs/<job_id>/logger/log.txt` 是否在写入。

### 9.6 大数据时“校验/Seurat 解析”长时间无响应

现版本前端已将 API 等待时间统一提升到 10 分钟（600 秒）。

如果数据很大仍感觉“卡住”，优先按下面检查：
1. 打开浏览器 F12 Network，确认对应请求是否还在 `pending`。
2. 查看后端终端是否仍在打印处理日志（尤其是 R 转换与 scanpy 读取阶段）。
3. 若最终返回超时提示，通常表示后端仍在算，可稍后重试，或先用更小数据验证流程。
4. Windows + Conda 场景建议保持 `cmd_conda`，并使用 `condabin\\conda.bat` 的完整路径。

---

## 10. 回滚说明

若要关闭新首页交互（自动进入流程页等），请回滚：
- `frontend/src/App.tsx`
- `frontend/src/styles/tokens.css`

若要关闭认证保护：
- 设置环境变量：`LABFLOW_AUTH_REQUIRED=false`
- 重启后端
