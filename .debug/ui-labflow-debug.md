# UI Debug: LabFlow 前端

## 设计决策记录

### [2026-02-11] 点击按钮后的即时反馈与运行日志「小电视」

- **用户需求**：每一步点按钮后要有即时反馈；训练等长任务要有实时监测，如进度条或把后端 shell 当「小电视」展示。
- **方案**：
  1. **全局「当前任务」反馈**：当存在进行中步骤（`busyStep`）或任务处于 queued/running 时，在系统状态下方展示一个固定条带，包含：
     - 无限循环动画的进度条（indeterminate），表示「正在处理」；
     - 文案：上传中… / 校验中… / 解析 Seurat 中… / 500×500 预处理中… / 提交训练中… / 训练运行中…。
  2. **运行日志「小电视」**：在「6) 任务轮询状态」中，当存在 job 时展示「运行日志」区块：
     - 轮询 `GET /api/jobs/{job_id}/log`，每 2.5 秒刷新；
     - 内容在深色、等宽、可滚动的 `<pre>` 中展示（max-height 220px），风格类似终端；
     - 任务结束后保留最后一次拉取的日志，与下方「查看日志」一致。
- **代码**：
  - `frontend/src/App.tsx`：新增 `liveJobLog` 状态；轮询 `getJobLog` 的 useEffect（job 为 queued/running 时）；`taskLabelMap`、`showTaskFeedback`、`taskLabel`；插入「当前任务」条带与「运行日志」面板。
  - `frontend/src/styles/tokens.css`：`.task-feedback`、`.task-progress`（含 `task-progress-shift` 动画）、`.task-label`、`.live-log-box`、`.live-log-pre`。
- **Checkfix**：`npm run lint`、`npm run build` 通过。
- **后续可做**：若后端为 validate/inspect/prepare 提供进度或流式日志，可再接上进度或实时日志。

### [2026-02-11] 校验错误就近展示 + R/conda 参数菜单化
- **用户需求**：报错应在用户视线处（按钮旁）弹出并给出当前系统推荐参数；R 相关参数（conda.bat、R 环境名）应由后端探测后在前端用菜单选择，避免手输。
- **实现**：① 校验失败时仅设置 validateError（不设 globalError），在「开始校验」按钮旁用 step-error 区块展示错误与 validateRecommendation（Rscript/Conda 时推荐 cmd_conda + 完整路径 + r-4.3）；表单项 onChange 时清除 validateError。② 后端新增 GET /api/runtime/conda-envs（backend/app/api/runtime.py）：Windows 下扫描 PATH 与常见路径得到 conda.bat 候选，执行 conda env list（--json 或解析文本）得到环境名列表；支持 query conda_bat 只返回该路径下的 envs。③ 前端：请求 getCondaEnvs() 填充 condaBatCandidates、condaEnvsList；conda.bat 用 select（候选 +「其他（手动输入）」）+ 可选 input；Conda R 环境名用 select（conda_envs）；切换 conda.bat 时重新请求 envs。④ 样式：.step-actions、.step-error、.step-error-msg、.step-error-recommend。
- **Checkfix**：ruff、frontend lint/build 通过。

### [2026-02-11] 参数名后「?」悬停 1 秒显示说明气泡（傻瓜化）
- **用户需求**：每个参数名称后面放一个小问号，鼠标悬停约 1 秒后弹出气泡，解释该参数在单细胞数据里是什么、干什么用，把用户当非专业用户「喂到嘴里」。
- **方案**：
  1. **ParamTooltip 组件**（`frontend/src/components/ParamTooltip.tsx`）：渲染内联「?」图标；`onMouseEnter` 设 1000ms 定时器，`onMouseLeave` 清除并隐藏；定时到后显示气泡（`role="tooltip"`），气泡内可继续悬停以保持显示；气泡深色背景、白字、小三角指向触发点，max-width 320px。
  2. **文案**：同一文件内导出 `PARAM_TOOLTIPS`（Record<string, string>），为「数据名称」「主数据文件」「SMILES CSV」「使用药物结构模式」「R 执行方式」「Rscript 命令」「conda.bat 路径」「Conda R 环境名」「分组字段（group_column）」「类群字段（cluster_column）」「随机种子（seed）」「待筛选 clusters（逗号分隔）」以及 gene_size、output_dim、batch_size、lr 撰写简短、单细胞/LabFlow 语境下的说明。
  3. **接入**：在 `App.tsx` 每个对应表单项的 label 文字后插入 `<ParamTooltip text={PARAM_TOOLTIPS["…"]} />`（或 `.gene_size` 等）。
- **样式**：`tokens.css` 新增 `.param-tooltip-wrap`、`.param-tooltip-trigger`（圆形灰底「?」、hover 时 accent 色）、`.param-tooltip-bubble`（absolute、深色、圆角、三角箭头）。
- **Checkfix**：`npm run build` 通过；无新增 lint 报错。

### [2026-02-11] Seurat 解析 metadata 值显示：添加细胞数统计（类似 R table()）
- **用户需求**：点击 metadata 字段（如 celltype）后，不仅要显示有哪些分类，还要显示每个分类有多少个细胞（类似 R 的 `table()` 函数），便于用户了解数据分布。
- **问题**：之前只返回唯一值列表，用户看到"共 0 个"，且无法知道每个分类的细胞数。
- **方案**：
  1. **后端**（`backend/app/services/seurat_inspector.py`）：使用 pandas Series 的 `value_counts()` 统计每个值的计数（类似 R `table()`），返回结构从 `dict[str, list[str]]` 改为 `dict[str, list[dict[str, Any]]]`，每个元素包含 `{"value": str, "count": int}`；按计数降序、值升序排序；跳过 NaN/None/空字符串。
  2. **前端类型**（`frontend/src/services/api.ts`）：`metadata_column_values` 类型从 `Record<string, string[]>` 改为 `Record<string, Array<{ value: string; count: number }>>`。
  3. **前端显示**（`frontend/src/App.tsx`）：在 metadata-values-list 中，每个 chip 显示 `{value} ({count})`，count 用 accent 色、加粗；hover 时 tooltip 显示完整信息。
- **样式**：`.value-chip` 改为 `inline-flex` 支持 gap，`.value-count` 用 accent 色、加粗、稍小字号。
- **Checkfix**：后端 `ruff check` 通过；前端 `npm run build` 通过。

### [2026-02-11] 修复 Conda R 环境名下拉选项不显示的问题
- **用户反馈**：Conda R 环境名在前端刷新不出来选项，只能手动输入。
- **问题诊断**：
  1. 初始加载时调用 `getCondaEnvs()` 无参数，后端可能无法获取环境列表（即使找到了 conda.bat candidates）
  2. 后端逻辑：无 `conda_bat` 参数时用第一个 candidate 获取环境列表，但可能执行失败返回空列表
  3. 手动输入 conda.bat 路径时，没有自动获取对应的环境列表
- **修复方案**：
  1. **初始加载优化**：获取到 `conda_bat_candidates` 后，如果 `conda_envs` 为空，用第一个 candidate 重新调用 `getCondaEnvs(firstBat)` 获取环境列表
  2. **手动输入增强**：在手动输入 conda.bat 路径的 `onChange` 中，如果输入以 "conda.bat" 结尾，自动调用 `getCondaEnvs()` 获取环境列表
- **代码**：`frontend/src/App.tsx` 第 87-105 行（初始加载 useEffect）和第 523-535 行（手动输入 input onChange）
- **Checkfix**：前端 `npm run build` 通过；无 lint 错误。

### [2026-02-11] ѵ����־��ѯ 404 �޸� + ǰ��ֹͣѵ������
- �û�������
  1. ѵ�����ܣ���ǰ�˳�����ʾ����ȡ��־ʧ�ܡ���
  2. ��Ҫ��ǰ��ֱ��ֹͣѵ������
- ����
  1. ����ѵ������ǰ�Ż�д `job.log_path`����ѯ `/api/jobs/{id}/log` ʱ���õ� 404��
  2. ���ȱ��ȡ������ API��ǰ��Ҳû�С�ֹͣ������ڡ�
- ʵ�֣�
  1. ��� `GET /api/jobs/{job_id}/log` ��Ϊ����־δ����Ҳ���� 200 + ���ַ�������
  2. ������� `POST /api/jobs/{job_id}/cancel`��queued ֱ�ӱ�� `canceled`��running ��� `cancel_requested` ��������ֹѵ�����̣�Windows �� `taskkill /T /F`����
  3. ���������ѵ������ǰԤд `log_path`�����ڲ���ȡ������ʱ��״̬����Ϊ `canceled`��
  4. ǰ��������ֹͣѵ�����񡱰�ť���� queued/running ��ʾ�������� cancel API����־��ѯʧ��ʱ���ٷ���д����ȡ��־ʧ�ܡ�ռλ�İ���
- �ĵ�ͬ����
  1. ���� `docs/api/jobs.md`������ cancel �� log ��Ϊ����
  2. ���� `docs/LabFlowǰ���û�����˵��.md`��������ֹͣѵ�����񡱲���˵������
- Checkfix��
  - ��ִ�в����
- Checkfix ������
  - `python -m black . --check`��ʧ�ܣ���ʾ 2 ���ļ����ʽ����`backend/app/api/runtime.py`��`backend/app/services/seurat_inspector.py`��
  - `npm run lint`��frontend����ͨ����
  - `npm run build`��frontend����ͨ����
  - ��� pytest����ǰ����ȱ�� pytest��`No module named pytest`����`uv run pytest` �ܱ���Ȩ�����ƣ��ܾ����ʣ���
- [2026-02-11] Checkfix ׷�ӣ�code-debugger��
  - �û�����ȷ�ϣ�`python -m black backend/app/services/seurat_inspector.py --check` ͨ����
  - ������Զ���ʽ��������`backend/app/api/runtime.py`��`backend/app/api/jobs.py`��ͳһ��ʽ�뻻�У���
  - ����ͨ����
    - `python -m black backend/app/api/runtime.py backend/app/services/seurat_inspector.py --check`
    - `ruff format --check backend/app backend/tests`
    - `ruff check backend/app backend/tests`
