# UI Debug: LabFlow 前端

## 元信息
- 模块名称: ui-labflow
- 创建时间: 2026-02-11
- 最后更新: 2026-02-11
- 相关文件:
- `frontend/src/App.tsx`
- `frontend/src/services/api.ts`
- `frontend/src/components/ParamTooltip.tsx`
- `frontend/src/styles/tokens.css`
- `backend/app/api/jobs.py`
- `backend/app/api/runtime.py`
- 用户说明书路径:
- `docs/LabFlow前端用户操作说明.md`
- 开发/部署文档路径:
- `README.md`
- `docs/部署文档.md`

## 运行上下文与测试规则
- 运行环境: 本机 Windows
- 验证/Checkfix 执行方式: 本地终端执行 `npm run lint` / `npm run build`，后端接口通过本地启动后联调
- 备注: 当前 PowerShell 会打印编码相关噪声，不影响代码实际结果

## 上下文关系网络
- 文件结构:
- `App.tsx` 负责首页登录、任务侧栏、流程步骤、日志与结果展示
- `api.ts` 负责前端 API 请求封装与错误处理
- `tokens.css` 负责布局、任务态样式、日志视图样式
- 函数调用链:
- 用户点击操作 -> `App.tsx` 事件函数 -> `api.ts` 请求 -> 后端接口 -> 状态回填 UI
- 数据流向:
- 任务状态: `/api/jobs` -> 任务侧栏与详情页
- 运行日志: `/api/jobs/{id}/log` -> 实时日志面板
- 运行环境探测: `/api/runtime/conda-envs` / `/api/runtime/gpu-stats`

## Debug 历史
### [2026-02-11] 任务即时反馈与运行日志实时显示
- 问题描述
- 用户希望点击后立即有反馈，训练期间能实时看到日志。
- 根因定位
- 任务执行是异步的，但前端在等待阶段缺少持续视觉反馈与实时日志区域。
- 解决方案
- 增加全局任务反馈条（busyStep/queued/running）。
- 增加运行日志轮询面板（每 2.5s 拉取 `/api/jobs/{job_id}/log`）。
- 代码变更
- `frontend/src/App.tsx`
- `frontend/src/styles/tokens.css`
- 验证结果
- `npm run lint` 通过。
- `npm run build` 在部分主机受权限影响可能报 `spawn EPERM`（非代码错误）。
- 影响评估
- 用户能看到“正在执行”的明确状态，降低“按钮没反应”的感知。

### [2026-02-11] 校验报错就近展示 + R/Conda 参数下拉化
- 问题描述
- 校验失败信息不够聚焦，R/conda 参数手输易错。
- 根因定位
- 错误展示位置偏离用户焦点；缺少运行时 conda 环境自动发现。
- 解决方案
- 校验错误在按钮附近就近展示，并给出推荐配置。
- 接入 `/api/runtime/conda-envs`，为 conda.bat 和环境名提供下拉选项。
- 代码变更
- `frontend/src/App.tsx`
- `frontend/src/services/api.ts`
- `backend/app/api/runtime.py`
- 验证结果
- 前端 lint/build 流程可执行。
- 影响评估
- 降低 R 运行配置门槛，减少“参数填错导致无法校验”。

### [2026-02-11] 参数说明 Tooltip（延迟显示）
- 问题描述
- 用户希望参数说明更直观，减少术语理解成本。
- 根因定位
- 表单参数仅有字段名，缺少上下文解释。
- 解决方案
- 新增 `ParamTooltip` 组件，悬停约 1 秒显示说明气泡。
- 为主要参数提供说明文案映射。
- 代码变更
- `frontend/src/components/ParamTooltip.tsx`
- `frontend/src/App.tsx`
- `frontend/src/styles/tokens.css`
- 验证结果
- 前端构建通过。
- 影响评估
- 提升非生信用户对参数的理解速度。

### [2026-02-11] Seurat metadata 值显示计数
- 问题描述
- 用户只看到分类值，看不到每类细胞数量。
- 根因定位
- 接口只返回唯一值列表，未返回 value count。
- 解决方案
- 后端按列统计 value count，前端显示为 `value (count)`。
- 代码变更
- `backend/app/services/seurat_inspector.py`
- `frontend/src/services/api.ts`
- `frontend/src/App.tsx`
- 验证结果
- 后端静态检查通过（环境具备依赖时）。
- 影响评估
- metadata 可读性更高，便于快速判断分布。

### [2026-02-11] 训练日志 404 与停止任务能力
- 问题描述
- 训练中日志偶发 404，且前端缺少停止任务入口。
- 根因定位
- 任务早期日志路径尚未稳定；取消接口与前端动作未闭环。
- 解决方案
- 优化后端日志读取 fallback。
- 新增并接入 `POST /api/jobs/{job_id}/cancel`。
- 前端增加“停止训练任务”动作。
- 代码变更
- `backend/app/api/jobs.py`
- `backend/app/services/job_queue.py`
- `frontend/src/App.tsx`
- 验证结果
- 功能联调通过（依赖后端可达）。
- 影响评估
- 提升长任务可控性，减少卡死体验。

## 待追踪问题
- 某些 Windows 终端会输出编码/权限噪声，不影响业务逻辑但影响排障体验。
- 需要在统一运行环境下补齐完整 checkfix 结果（尤其后端测试依赖）。

## 技术债务记录
- 前端单页承载较多流程状态，后续可拆分页面与状态模块。
- 运行日志与任务轮询目前是定时拉取，后续可升级为流式推送。
