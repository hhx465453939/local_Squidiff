# UI Debug: LabFlow 鍓嶇

## 璁捐鍐崇瓥璁板綍

### [2026-02-11] 鐐瑰嚮鎸夐挳鍚庣殑鍗虫椂鍙嶉涓庤繍琛屾棩蹇椼€屽皬鐢佃銆?

- **鐢ㄦ埛闇€姹?*锛氭瘡涓€姝ョ偣鎸夐挳鍚庤鏈夊嵆鏃跺弽棣堬紱璁粌绛夐暱浠诲姟瑕佹湁瀹炴椂鐩戞祴锛屽杩涘害鏉℃垨鎶婂悗绔?shell 褰撱€屽皬鐢佃銆嶅睍绀恒€?
- **鏂规**锛?
  1. **鍏ㄥ眬銆屽綋鍓嶄换鍔°€嶅弽棣?*锛氬綋瀛樺湪杩涜涓楠わ紙`busyStep`锛夋垨浠诲姟澶勪簬 queued/running 鏃讹紝鍦ㄧ郴缁熺姸鎬佷笅鏂瑰睍绀轰竴涓浐瀹氭潯甯︼紝鍖呭惈锛?
     - 鏃犻檺寰幆鍔ㄧ敾鐨勮繘搴︽潯锛坕ndeterminate锛夛紝琛ㄧず銆屾鍦ㄥ鐞嗐€嶏紱
     - 鏂囨锛氫笂浼犱腑鈥?/ 鏍￠獙涓€?/ 瑙ｆ瀽 Seurat 涓€?/ 500脳500 棰勫鐞嗕腑鈥?/ 鎻愪氦璁粌涓€?/ 璁粌杩愯涓€︺€?
  2. **杩愯鏃ュ織銆屽皬鐢佃銆?*锛氬湪銆?) 浠诲姟杞鐘舵€併€嶄腑锛屽綋瀛樺湪 job 鏃跺睍绀恒€岃繍琛屾棩蹇椼€嶅尯鍧楋細
     - 杞 `GET /api/jobs/{job_id}/log`锛屾瘡 2.5 绉掑埛鏂帮紱
     - 鍐呭鍦ㄦ繁鑹层€佺瓑瀹姐€佸彲婊氬姩鐨?`<pre>` 涓睍绀猴紙max-height 220px锛夛紝椋庢牸绫讳技缁堢锛?
     - 浠诲姟缁撴潫鍚庝繚鐣欐渶鍚庝竴娆℃媺鍙栫殑鏃ュ織锛屼笌涓嬫柟銆屾煡鐪嬫棩蹇椼€嶄竴鑷淬€?
- **浠ｇ爜**锛?
  - `frontend/src/App.tsx`锛氭柊澧?`liveJobLog` 鐘舵€侊紱杞 `getJobLog` 鐨?useEffect锛坖ob 涓?queued/running 鏃讹級锛沗taskLabelMap`銆乣showTaskFeedback`銆乣taskLabel`锛涙彃鍏ャ€屽綋鍓嶄换鍔°€嶆潯甯︿笌銆岃繍琛屾棩蹇椼€嶉潰鏉裤€?
  - `frontend/src/styles/tokens.css`锛歚.task-feedback`銆乣.task-progress`锛堝惈 `task-progress-shift` 鍔ㄧ敾锛夈€乣.task-label`銆乣.live-log-box`銆乣.live-log-pre`銆?
- **Checkfix**锛歚npm run lint`銆乣npm run build` 閫氳繃銆?
- **鍚庣画鍙仛**锛氳嫢鍚庣涓?validate/inspect/prepare 鎻愪緵杩涘害鎴栨祦寮忔棩蹇楋紝鍙啀鎺ヤ笂杩涘害鎴栧疄鏃舵棩蹇椼€?

### [2026-02-11] 鏍￠獙閿欒灏辫繎灞曠ず + R/conda 鍙傛暟鑿滃崟鍖?
- **鐢ㄦ埛闇€姹?*锛氭姤閿欏簲鍦ㄧ敤鎴疯绾垮锛堟寜閽梺锛夊脊鍑哄苟缁欏嚭褰撳墠绯荤粺鎺ㄨ崘鍙傛暟锛汻 鐩稿叧鍙傛暟锛坈onda.bat銆丷 鐜鍚嶏級搴旂敱鍚庣鎺㈡祴鍚庡湪鍓嶇鐢ㄨ彍鍗曢€夋嫨锛岄伩鍏嶆墜杈撱€?
- **瀹炵幇**锛氣憼 鏍￠獙澶辫触鏃朵粎璁剧疆 validateError锛堜笉璁?globalError锛夛紝鍦ㄣ€屽紑濮嬫牎楠屻€嶆寜閽梺鐢?step-error 鍖哄潡灞曠ず閿欒涓?validateRecommendation锛圧script/Conda 鏃舵帹鑽?cmd_conda + 瀹屾暣璺緞 + r-4.3锛夛紱琛ㄥ崟椤?onChange 鏃舵竻闄?validateError銆傗憽 鍚庣鏂板 GET /api/runtime/conda-envs锛坆ackend/app/api/runtime.py锛夛細Windows 涓嬫壂鎻?PATH 涓庡父瑙佽矾寰勫緱鍒?conda.bat 鍊欓€夛紝鎵ц conda env list锛?-json 鎴栬В鏋愭枃鏈級寰楀埌鐜鍚嶅垪琛紱鏀寔 query conda_bat 鍙繑鍥炶璺緞涓嬬殑 envs銆傗憿 鍓嶇锛氳姹?getCondaEnvs() 濉厖 condaBatCandidates銆乧ondaEnvsList锛沜onda.bat 鐢?select锛堝€欓€?+銆屽叾浠栵紙鎵嬪姩杈撳叆锛夈€嶏級+ 鍙€?input锛汣onda R 鐜鍚嶇敤 select锛坈onda_envs锛夛紱鍒囨崲 conda.bat 鏃堕噸鏂拌姹?envs銆傗懀 鏍峰紡锛?step-actions銆?step-error銆?step-error-msg銆?step-error-recommend銆?
- **Checkfix**锛歳uff銆乫rontend lint/build 閫氳繃銆?

### [2026-02-11] 鍙傛暟鍚嶅悗銆?銆嶆偓鍋?1 绉掓樉绀鸿鏄庢皵娉★紙鍌荤摐鍖栵級
- **鐢ㄦ埛闇€姹?*锛氭瘡涓弬鏁板悕绉板悗闈㈡斁涓€涓皬闂彿锛岄紶鏍囨偓鍋滅害 1 绉掑悗寮瑰嚭姘旀场锛岃В閲婅鍙傛暟鍦ㄥ崟缁嗚優鏁版嵁閲屾槸浠€涔堛€佸共浠€涔堢敤锛屾妸鐢ㄦ埛褰撻潪涓撲笟鐢ㄦ埛銆屽杺鍒板槾閲屻€嶃€?
- **鏂规**锛?
  1. **ParamTooltip 缁勪欢**锛坄frontend/src/components/ParamTooltip.tsx`锛夛細娓叉煋鍐呰仈銆?銆嶅浘鏍囷紱`onMouseEnter` 璁?1000ms 瀹氭椂鍣紝`onMouseLeave` 娓呴櫎骞堕殣钘忥紱瀹氭椂鍒板悗鏄剧ず姘旀场锛坄role="tooltip"`锛夛紝姘旀场鍐呭彲缁х画鎮仠浠ヤ繚鎸佹樉绀猴紱姘旀场娣辫壊鑳屾櫙銆佺櫧瀛椼€佸皬涓夎鎸囧悜瑙﹀彂鐐癸紝max-width 320px銆?
  2. **鏂囨**锛氬悓涓€鏂囦欢鍐呭鍑?`PARAM_TOOLTIPS`锛圧ecord<string, string>锛夛紝涓恒€屾暟鎹悕绉般€嶃€屼富鏁版嵁鏂囦欢銆嶃€孲MILES CSV銆嶃€屼娇鐢ㄨ嵂鐗╃粨鏋勬ā寮忋€嶃€孯 鎵ц鏂瑰紡銆嶃€孯script 鍛戒护銆嶃€宑onda.bat 璺緞銆嶃€孋onda R 鐜鍚嶃€嶃€屽垎缁勫瓧娈碉紙group_column锛夈€嶃€岀被缇ゅ瓧娈碉紙cluster_column锛夈€嶃€岄殢鏈虹瀛愶紙seed锛夈€嶃€屽緟绛涢€?clusters锛堥€楀彿鍒嗛殧锛夈€嶄互鍙?gene_size銆乷utput_dim銆乥atch_size銆乴r 鎾板啓绠€鐭€佸崟缁嗚優/LabFlow 璇涓嬬殑璇存槑銆?
  3. **鎺ュ叆**锛氬湪 `App.tsx` 姣忎釜瀵瑰簲琛ㄥ崟椤圭殑 label 鏂囧瓧鍚庢彃鍏?`<ParamTooltip text={PARAM_TOOLTIPS["鈥?]} />`锛堟垨 `.gene_size` 绛夛級銆?
- **鏍峰紡**锛歚tokens.css` 鏂板 `.param-tooltip-wrap`銆乣.param-tooltip-trigger`锛堝渾褰㈢伆搴曘€?銆嶃€乭over 鏃?accent 鑹诧級銆乣.param-tooltip-bubble`锛坅bsolute銆佹繁鑹层€佸渾瑙掋€佷笁瑙掔澶达級銆?
- **Checkfix**锛歚npm run build` 閫氳繃锛涙棤鏂板 lint 鎶ラ敊銆?

### [2026-02-11] Seurat 瑙ｆ瀽 metadata 鍊兼樉绀猴細娣诲姞缁嗚優鏁扮粺璁★紙绫讳技 R table()锛?
- **鐢ㄦ埛闇€姹?*锛氱偣鍑?metadata 瀛楁锛堝 celltype锛夊悗锛屼笉浠呰鏄剧ず鏈夊摢浜涘垎绫伙紝杩樿鏄剧ず姣忎釜鍒嗙被鏈夊灏戜釜缁嗚優锛堢被浼?R 鐨?`table()` 鍑芥暟锛夛紝渚夸簬鐢ㄦ埛浜嗚В鏁版嵁鍒嗗竷銆?
- **闂**锛氫箣鍓嶅彧杩斿洖鍞竴鍊煎垪琛紝鐢ㄦ埛鐪嬪埌"鍏?0 涓?锛屼笖鏃犳硶鐭ラ亾姣忎釜鍒嗙被鐨勭粏鑳炴暟銆?
- **鏂规**锛?
  1. **鍚庣**锛坄backend/app/services/seurat_inspector.py`锛夛細浣跨敤 pandas Series 鐨?`value_counts()` 缁熻姣忎釜鍊肩殑璁℃暟锛堢被浼?R `table()`锛夛紝杩斿洖缁撴瀯浠?`dict[str, list[str]]` 鏀逛负 `dict[str, list[dict[str, Any]]]`锛屾瘡涓厓绱犲寘鍚?`{"value": str, "count": int}`锛涙寜璁℃暟闄嶅簭銆佸€煎崌搴忔帓搴忥紱璺宠繃 NaN/None/绌哄瓧绗︿覆銆?
  2. **鍓嶇绫诲瀷**锛坄frontend/src/services/api.ts`锛夛細`metadata_column_values` 绫诲瀷浠?`Record<string, string[]>` 鏀逛负 `Record<string, Array<{ value: string; count: number }>>`銆?
  3. **鍓嶇鏄剧ず**锛坄frontend/src/App.tsx`锛夛細鍦?metadata-values-list 涓紝姣忎釜 chip 鏄剧ず `{value} ({count})`锛宑ount 鐢?accent 鑹层€佸姞绮楋紱hover 鏃?tooltip 鏄剧ず瀹屾暣淇℃伅銆?
- **鏍峰紡**锛歚.value-chip` 鏀逛负 `inline-flex` 鏀寔 gap锛宍.value-count` 鐢?accent 鑹层€佸姞绮椼€佺◢灏忓瓧鍙枫€?
- **Checkfix**锛氬悗绔?`ruff check` 閫氳繃锛涘墠绔?`npm run build` 閫氳繃銆?

### [2026-02-11] 淇 Conda R 鐜鍚嶄笅鎷夐€夐」涓嶆樉绀虹殑闂
- **鐢ㄦ埛鍙嶉**锛欳onda R 鐜鍚嶅湪鍓嶇鍒锋柊涓嶅嚭鏉ラ€夐」锛屽彧鑳芥墜鍔ㄨ緭鍏ャ€?
- **闂璇婃柇**锛?
  1. 鍒濆鍔犺浇鏃惰皟鐢?`getCondaEnvs()` 鏃犲弬鏁帮紝鍚庣鍙兘鏃犳硶鑾峰彇鐜鍒楄〃锛堝嵆浣挎壘鍒颁簡 conda.bat candidates锛?
  2. 鍚庣閫昏緫锛氭棤 `conda_bat` 鍙傛暟鏃剁敤绗竴涓?candidate 鑾峰彇鐜鍒楄〃锛屼絾鍙兘鎵ц澶辫触杩斿洖绌哄垪琛?
  3. 鎵嬪姩杈撳叆 conda.bat 璺緞鏃讹紝娌℃湁鑷姩鑾峰彇瀵瑰簲鐨勭幆澧冨垪琛?
- **淇鏂规**锛?
  1. **鍒濆鍔犺浇浼樺寲**锛氳幏鍙栧埌 `conda_bat_candidates` 鍚庯紝濡傛灉 `conda_envs` 涓虹┖锛岀敤绗竴涓?candidate 閲嶆柊璋冪敤 `getCondaEnvs(firstBat)` 鑾峰彇鐜鍒楄〃
  2. **鎵嬪姩杈撳叆澧炲己**锛氬湪鎵嬪姩杈撳叆 conda.bat 璺緞鐨?`onChange` 涓紝濡傛灉杈撳叆浠?"conda.bat" 缁撳熬锛岃嚜鍔ㄨ皟鐢?`getCondaEnvs()` 鑾峰彇鐜鍒楄〃
- **浠ｇ爜**锛歚frontend/src/App.tsx` 绗?87-105 琛岋紙鍒濆鍔犺浇 useEffect锛夊拰绗?523-535 琛岋紙鎵嬪姩杈撳叆 input onChange锛?
- **Checkfix**锛氬墠绔?`npm run build` 閫氳繃锛涙棤 lint 閿欒銆?

### [2026-02-11] 训练日志轮询 404 修复 + 前端停止训练任务
- 用户反馈：
  1. 训练在跑，但前端持续显示“获取日志失败”。
  2. 需要在前端直接停止训练任务。
- 根因：
  1. 任务训练结束前才回写 `job.log_path`，轮询 `/api/jobs/{id}/log` 时常拿到 404。
  2. 后端缺少取消任务 API，前端也没有“停止任务”入口。
- 实现：
  1. 后端 `GET /api/jobs/{job_id}/log` 改为“日志未就绪也返回 200 + 空字符串”。
  2. 后端新增 `POST /api/jobs/{job_id}/cancel`：queued 直接标记 `canceled`；running 标记 `cancel_requested` 并尝试终止训练进程（Windows 用 `taskkill /T /F`）。
  3. 任务队列在训练启动前预写 `log_path`，并在捕获到取消请求时将状态收敛为 `canceled`。
  4. 前端新增“停止训练任务”按钮（仅 queued/running 显示），接入 cancel API；日志轮询失败时不再反复写“获取日志失败”占位文案。
- 文档同步：
  1. 新增 `docs/api/jobs.md`（包含 cancel 与 log 行为）。
  2. 更新 `docs/LabFlow前端用户操作说明.md`（新增“停止训练任务”操作说明）。
- Checkfix：
  - 待执行并回填。
- Checkfix 结果回填：
  - `python -m black . --check`：失败，提示 2 个文件需格式化：`backend/app/api/runtime.py`、`backend/app/services/seurat_inspector.py`。
  - `npm run lint`（frontend）：通过。
  - `npm run build`（frontend）：通过。
  - 后端 pytest：当前环境缺少 pytest（`No module named pytest`）；`uv run pytest` 受本机权限限制（拒绝访问）。
- [2026-02-11] Checkfix 追加（code-debugger）
  - 用户本机确认：`python -m black backend/app/services/seurat_inspector.py --check` 通过。
  - 已完成自动格式化收敛：`backend/app/api/runtime.py`、`backend/app/api/jobs.py`（统一格式与换行）。
  - 复检通过：
    - `python -m black backend/app/api/runtime.py backend/app/services/seurat_inspector.py --check`
    - `ruff format --check backend/app backend/tests`
    - `ruff check backend/app backend/tests`

### [2026-02-11] 动画首页 + 开始分析入口（科研视觉升维）
- 艺术指导：
  - Mood: 理性、精密、可观测。
  - Metaphor: “实验台上的光学扫描屏”，先展示系统能力，再进入分析流程。
- 视觉审计与策略：
  1. 空间：原页面打开即进入多面板，信息密度高；新增 Landing 首屏，分层展示能力点与入口动作。
  2. 张力：加入扫描高光、网格纹理、卡片呼吸动画，增强科研仪表感。
  3. 质感：多层渐变 + 半透明卡片 + 轻阴影，避免纯平界面。
  4. 微交互：入口按钮 hover 提升、流程页进入时上浮过渡。
- 实施记录：
  - `frontend/src/App.tsx`：新增 `hasEntered` 入口状态；加入动画首页 DOM；点击“开始分析”后渲染原有 1~7 流程。
  - `frontend/src/styles/tokens.css`：新增 landing 相关 token 与样式（网格层、扫描动画、卡片呼吸、入口按钮动效、移动端适配）。
- 用户说明书更新：
  - `docs/LabFlow前端用户操作说明.md` 新增“动画首页入口（开始分析）”章节（目的、操作、预期、常见问题、回滚）。
- 部署文档联动检查：
  - 已检查 `docs/部署文档.md`，本次仅前端视觉与入口交互变化，不涉及部署命令与环境变量，无需改动。
- Checkfix：
  - `npm run lint`（frontend）通过。
  - `npm run build`（frontend）通过。

### [2026-02-11] 首页注册/登录 + 用户说明书入口（内网轻量认证）
- 目标：在动画首页加入可直接使用的账号入口，登录后进入分析流程，并把用户说明书按钮放在首页。
- 后端实现（API-First）：
  - 新增 `backend/app/services/auth_service.py`：SQLite 本地库 + PBKDF2-SHA256 密码哈希 + session token。
  - 新增 `backend/app/api/auth.py`：`/api/auth/register`、`/api/auth/login`、`/api/auth/me`、`/api/auth/logout`、`/api/auth/user-guide`。
  - 新增 `backend/app/auth.py`：Bearer 解析与认证依赖。
  - `backend/app/main.py` 路由挂载：业务 API 统一走 `require_auth`（可由 `LABFLOW_AUTH_REQUIRED` 控制）。
  - `backend/app/core/config.py` 增加认证配置：`LABFLOW_AUTH_REQUIRED`、`LABFLOW_AUTH_DB_PATH`、`LABFLOW_AUTH_SESSION_TTL_HOURS`。
- 前端实现（首页入口与交互）：
  - `frontend/src/App.tsx`：首页新增注册/登录表单、登录态显示、退出登录、未登录禁用“开始分析”。
  - 首页“用户说明书”按钮改为后端在线文档接口：`/api/auth/user-guide`。
  - `frontend/src/services/api.ts`：新增 auth API 封装与 token 本地存取；请求自动携带 Bearer Token。
  - `frontend/src/styles/tokens.css`：新增 auth panel、tabs、按钮样式。
- 文档同步：
  - 新增 `docs/api/auth.md`。
  - 更新 `docs/LabFlow前端用户操作说明.md`（新增首页注册/登录步骤与规则）。
  - 更新 `docs/部署文档.md`（新增认证相关环境变量）。
  - 更新 `README.md`（新增登录认证说明入口）。
- 用户可理解性检查结论：
  - 结论：基本可让新手上手。说明书已覆盖“先登录再开始分析”的路径、账号规则、常见报错。
  - 仍建议（后续可做）：补一张“注册->登录->开始分析”流程截图，进一步降低首次认知成本。
- Checkfix：
  - `ruff check backend/app backend/tests` 通过。
  - `ruff format --check backend/app backend/tests` 通过。
  - `npm run lint`（frontend）通过。
  - `npm run build`（frontend）通过。
  - `python -m pytest -q backend/tests/test_auth_api.py` 未执行（当前环境缺少 pytest）。
  - 运行时 smoke（设置 `LABFLOW_AUTH_DB_PATH=.debug/auth_smoke.db`）通过：
    - `POST /api/auth/register` -> 200
    - `GET /api/auth/me` -> 200
    - `POST /api/auth/logout` -> 200
    - logout 后 `GET /api/auth/me` -> 401
    - `GET /api/auth/user-guide` -> 200
    - 未登录访问 `GET /api/jobs` -> 401

### [2026-02-11 16:31 +08:00] 任务轮询体验增强：GPU 实时看板 + 登录直达流程 + 日志早显
- 问题描述
  - 用户希望轮询阶段有更直观反馈（动画或显卡状态），并且日志不要等“停止任务”后才出现。
  - 登录后仍需再点“开始分析”，流程不够顺滑。
  - 数据提交页缺少用户说明书入口。
- 根因定位
  - 前端未接入 GPU 轮询数据源。
  - 登录成功只设置 token/user，未切换 `hasEntered=true`。
  - 任务日志接口优先读取 `train.log`，训练早期该文件可能为空，导致前端显示“等待日志”。
- 解决方案
  - 后端新增 `GET /api/runtime/gpu-stats`（`nvidia-smi` 解析）。
  - 前端在任务 `running` 时轮询 GPU 指标并展示 `GPU Runtime` 卡片（含脉冲动画和资源条）。
  - 登录成功后自动 `setHasEntered(true)` 直达分析页。
  - 在 `1) 上传数据` 面板新增 `User Guide` 按钮。
  - `GET /api/jobs/{job_id}/log` 增加 fallback：当 `job.log_path` 内容为空时，回退读取 `backend/artifacts/jobs/{job_id}/logger/log.txt`。
  - 停止训练按钮仅在 `running` 显示，不在 `queued` 显示。
- 代码变更（文件）
  - `backend/app/api/runtime.py`
  - `backend/app/api/jobs.py`
  - `frontend/src/services/api.ts`
  - `frontend/src/App.tsx`
  - `frontend/src/styles/tokens.css`
  - `docs/api/runtime.md`（新增）
  - `docs/api/auth.md`
  - `docs/api/jobs.md`
  - `docs/LabFlow前端用户操作说明.md`
- Checkfix 结果
  - `ruff format --check backend/app backend/tests` -> 通过（先自动格式化 `backend/app/api/auth.py`）
  - `ruff check backend/app backend/tests` -> 通过
  - `npm run lint`（frontend）-> 通过
  - `npm run build`（frontend）-> 通过
- 部署/文档联动检查
  - 已检查部署文档：本次不新增部署命令，仅新增运行时 API 与前端展示，无需修改 `docs/部署文档.md`。

### [2026-02-11 16:45 +08:00] 淇锛?api/auth/user-guide 鎵撲笉寮€ + 鏂囨。涔辩爜 + 鍓嶇杩炴帴鎷掔粷
- 闂鐜拌薄
  - 鍓嶇鎺у埗鍙板ぇ閲?`ERR_CONNECTION_REFUSED`锛坄/api/health`銆乣/api/auth/me`銆乣/api/auth/register`锛夈€?  - `http://localhost:8000/api/auth/user-guide` 鏃犳硶璁块棶銆?  - `docs/LabFlow鍓嶇鐢ㄦ埛鎿嶄綔璇存槑.md` 鍚庡崐娈典贡鐮併€?- 鏍瑰洜
  1. 鍚庣鍚姩澶辫触锛歚backend/app/api/auth.py` 涓?`user_guide` 杩斿洖绫诲瀷鏍囨敞涓?`FileResponse | HTMLResponse`锛孎astAPI 鍦ㄨ矾鐢卞缓妯℃椂鎶ラ敊骞堕€€鍑恒€?  2. 鐢ㄦ埛璇存槑涔︽枃浠剁紪鐮佸凡鎹熷潖锛屽鑷?`/api/auth/user-guide` 杩斿洖 `User guide encoding is not supported`銆?- 淇
  1. `backend/app/api/auth.py`
     - `@router.get("/user-guide", response_model=None)`
     - 杩斿洖绫诲瀷鏀逛负 `Response`锛岄伩鍏?FastAPI 灏嗗搷搴旂被 Union 褰撴垚 Pydantic 瀛楁銆?  2. 褰诲簳閲嶅缓 `docs/LabFlow鍓嶇鐢ㄦ埛鎿嶄綔璇存槑.md`
     - 浠?UTF-8 閲嶅啓鏁翠唤鏂囨。锛岃鐩栫櫥褰曘€佷笂浼犮€佹牎楠屻€佽缁冦€丟PU 鐪嬫澘銆佸父瑙佹晠闅滀笌鍥炴粴銆?- 楠岃瘉
  - 鏈湴鍚姩鍚庣骞惰皟鐢細
    - `GET /api/health` -> 200
    - `GET /api/auth/user-guide` -> 200锛圚TML锛?    - `GET /api/auth/user-guide?raw=true` -> 200锛圡arkdown锛?  - Checkfix锛?    - `ruff format --check backend/app backend/tests` 閫氳繃
    - `ruff check backend/app backend/tests` 閫氳繃
    - `npm run lint` 閫氳繃
    - `npm run build` 閫氳繃
- 鐢ㄦ埛渚ф搷浣滄彁绀?  - 鍑虹幇 `ERR_CONNECTION_REFUSED` 鏃跺厛纭鍚庣宸茶繍琛屽湪 `8000` 绔彛銆?

### [2026-02-11 17:20 +08:00] README 与 docs 编码乱码修复
- 问题描述
  - 仓库首页 README 出现乱码，要求同步排查 docs 文档。
- 根因定位
  - 文档编码在不同环境下被错误识别（UTF-8 无 BOM 在部分环境下被当作其他编码解析）。
- 处理动作
  - 统一将 `README.md` 与 `docs/**/*.md` 重写为 UTF-8 BOM。
  - 全量复检 UTF-8 严格解码通过。
  - 清理遗留备份文件：`docs/部署文档.md.__bak__`。
- 验证结果
  - `README.md`、`docs/**/*.md` 均可 UTF-8 严格解码。
  - 所有目标文件 BOM 检查通过。
- 影响评估
  - 仅编码规范化，不涉及功能逻辑变更。
