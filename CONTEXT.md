# GeoRBF

GeoRBF 描述以地质观测约束隐式标量场的领域语言，使使用者不必接触 RBF 装配与求解器细节。

## Language

**GeoRBF 使用者**:
掌握地质建模输入、并将 GeoRBF 嵌入 Rust 应用的开发者；不要求其具备 RBF 装配或数值求解器知识。
_Avoid_: RBF 专家、求解器使用者

**隐式地质标量场**:
由界面、方向、层位关系及不等式等地质观测共同约束，并可在空间位置上求取场值与梯度的标量场。
_Avoid_: 曲面模型、可视化模型

**层位（Horizon）**:
地层学语义中的一组空间位置，它们共享同一个未知层位场值；该共享场值属于问题语义。
_Avoid_: 已知等值面、参考点组、Interface

**地层场值方向（Stratigraphic Field Direction）**:
对一个地层标量场显式声明场值朝年轻方向增大或朝古老方向增大，从而把地层年龄关系映射为共享场值顺序。
_Avoid_: 默认场值方向、Z 轴方向、坐标手性、输入顺序

**地层年龄关系（Stratigraphic Age Relation）**:
使用 `YoungerThan` 或 `OlderThan` 陈述两个层位的相对地质年龄，并通过地层场值方向解释其数值顺序。
_Avoid_: Above、Below、空间高低关系

**最小场值间隔（Minimum Field Separation）**:
为严格地层年龄关系或命名关系分组显式规定的最小正共享场值差；它参与确定标量场尺度，单位为场值单位。
_Avoid_: 物理厚度、距离、数值 epsilon、求解容差

**场值间隔区间（Field Separation Interval）**:
以显式 reference 和 target 为有序角色，对两个共享等值集的带符号场值差规定有限双侧区间。
_Avoid_: Minimum Field Separation、地层年龄关系、物理厚度、无方向的点对

**物理厚度（Physical Thickness）**:
以长度单位表示、沿明确测量路径或规则定义的两个地质界面之间的几何距离；它不同于标量场中的场值间隔。
_Avoid_: Minimum Field Separation、Field Separation Interval、假设梯度模长后的场值换算

**共享等值集（Shared Level Set）**:
以稳定组身份和显式潜在共享场值关联一组空间位置的通用数学概念，可承载层位及其他等值集语义。
_Avoid_: 数学核心中的 Horizon、Interface、锚点组

**无信息共享等值集（Uninformative Shared Level Set）**:
只有一个成员、且潜在共享场值没有被任何其他关系引用，因而不向拟合问题提供约束信息的共享等值集。
_Avoid_: 被层序、间隔或 gauge 引用的单成员组

**场值层级关系（Field Level Order）**:
直接陈述两个共享等值集的未知场值大小关系，不附加地层年龄或空间高低含义。
_Avoid_: YoungerThan、OlderThan、Above、Below

**点—等值集关系（Point to Level Set Relation）**:
陈述一个指定位置相对某个显式共享等值集处于场值增加侧或场值减小侧的有限采样关系。
_Avoid_: off-interface inequality、区域侧保证、隐式参考界面

**最小场值偏移（Minimum Field Offset）**:
为严格点—等值集关系显式规定的最小正场值差，单位为场值单位。
_Avoid_: Minimum Field Separation、距离、物理厚度、数值 epsilon、求解容差

**场值观测（Field Value Observation）**:
在一个空间位置观测到已知绝对标量值的观测，与共享未知场值的层位或共享等值集不同。
_Avoid_: 已知层位、已知 isovalue 的 Horizon

**场值边界（Field Value Bound）**:
在一个指定位置对绝对场值规定显式有限的下界、上界或双侧区间；它会参与固定标量场的加性 gauge。
_Avoid_: Point to Level Set Relation、连续区域范围、NaN 或无穷哨兵

**加性场值 Gauge（Additive Field Gauge）**:
当全部领域关系对全局常数偏移不变时，显式选择一个绝对场值代表的约定；它不是测量观测，也不得改变梯度、场值差或等值集几何。
_Avoid_: 自动 anchor、首个层位为零、Field Value Observation、solver 最小范数

**未识别场模式（Unidentified Field Mode）**:
除全局常数偏移外，无法由当前领域观测确定、并会改变梯度或几何等可观察量的场变化模式。
_Avoid_: Polynomial Gauge、隐藏正则化、自动最小范数代表

**梯度观测（Gradient Observation）**:
在一个空间位置观测到完整梯度向量的观测；该向量同时包含方向、极性和梯度模长，单位为场值单位除以长度单位。
_Avoid_: Normal、方向观测、另附尺度的法向

**梯度残差（Gradient Residual）**:
预测梯度与完整梯度观测之差形成的物理坐标向量残差，以欧氏范数或显式协方差整体解释。
_Avoid_: 逐坐标 L1、方向与模长的隐式拆分

**法向方向观测（Normal Direction Observation）**:
在一个正则等值面点观测到非零梯度的几何法向方向，分为有向法向与无向法向轴；输入向量的原始模长没有领域意义。
_Avoid_: Gradient Observation、梯度模长、置信度或权重

**有向法向（Directed Normal）**:
包含极性的单位法向方向；它表示梯度为该方向的正倍数，因此 \(n\) 与 \(-n\) 不等价，零梯度不满足该观测。
_Avoid_: Axial Normal、无符号法向

**有向法向锥（Directed Normal Cone）**:
在输入坐标框架的欧氏度量中，把非零梯度限制在有向法向周围指定半角内的方向关系。
_Avoid_: Axial Normal 双锥、各向异性度量中的角度、分量 box

**法向锥违反量（Cone Violation）**:
软化有向法向锥边界的独立非负标量违反量，单位为场值单位除以长度单位。
_Avoid_: 角度误差、Minimum Normal Slope slack、共同法向 slack

**法向方向残差（Normal Direction Residual）**:
梯度在有向法向切平面上的旋转不变投影残差；它只衡量方向偏离，不衡量最小斜率是否满足。
_Avoid_: 任意切平面基的逐分量残差、完整法向满足状态

**角度（Angle）**:
以 `degree` 或 `radian` 显式标注单位、并在领域边界转换为规范弧度值的角度量。
_Avoid_: 裸标量、默认角度单位、根据数值猜测单位

**最小法向斜率（Minimum Normal Slope）**:
独立规定法向方向位置或命名观测分组沿所选法向的最小正梯度尺度，单位为场值单位除以长度单位。
_Avoid_: 法向量模长、权重、置信度、数值 epsilon、求解容差

**无向法向轴（Axial Normal）**:
不包含极性的法向方向等价类；它表示梯度沿该轴任一方向且不低于适用的最小法向斜率，因此 \(n\) 与 \(-n\) 等价。
_Avoid_: Directed Normal、带符号法向

**极性解析（Polarity Resolution）**:
把无向法向轴显式转换为有向法向的可追踪决定，记录所选符号及其来源。
_Avoid_: 启发式选符号、静默定向、自动翻转

**法向（Normal）**:
面向地质使用者对有向法向和无向法向轴的统称，不表示完整梯度观测。
_Avoid_: 底层含混观测类型、Gradient Observation

**切向方向观测（Tangent Direction Observation）**:
在一个空间位置观测到无极性的单位方向，并陈述场沿该方向的一阶方向导数为零；该观测允许零梯度，不能单独保证正则等值面。
_Avoid_: 非零梯度保证、带模长的切向、隐式法向选择

**方向导数区间（Directional Derivative Interval）**:
把场沿指定单位方向的一阶方向导数限制在一个有限区间内，区间单位为场值单位除以长度单位。
_Avoid_: 切向角度容差、假定梯度模长的角度换算、数值容差

**地质界面（Interface）**:
地质领域中的边界或接触面分类，包括层位、断层面和侵入接触面；分类本身不规定数值约束语义。
_Avoid_: 通用数值约束、Shared Level Set 的同义词

**输入坐标框架**:
调用方为一个问题显式提供的有序正交三维笛卡尔基、手性和公共长度单位；所有空间输入及输出都在该框架中解释。
_Avoid_: CRS、经纬度坐标、默认 Z 向上、内部归一化坐标、非均匀单位缩放

**坐标框架转换（Coordinate Frame Transform）**:
输入坐标框架之间由平移、正交旋转或反射及统一正尺度组成的显式可逆转换，所有方向、导数量与不确定性随之协变。
_Avoid_: 非均匀缩放、shear、静默单位转换、内部归一化结果

**全局各向异性度量（Global Anisotropy Metric）**:
一个问题内统一使用的无量纲、对称正定且行列式为一的度量，只表达核距离的方向与轴比，不改变物理坐标中的角度语义。
_Avoid_: 核长度尺度、局部 metric field、从法向自动估计、自动修正的半正定矩阵

**求解诊断**:
每次求解尝试产生的结构化结果，用于判别输入无效、欠约束、不可行或数值失败，并量化成功解的残差与约束违反。
_Avoid_: 日志、错误字符串、求解器原始输出

**求解尝试终止状态（Solve Attempt Termination）**:
一次具体后端调用停止的方式，包括成功候选、降低精度候选、不可行候选、达到限制或数值错误；它本身不证明问题性质。
_Avoid_: 问题诊断、已求解模型、后端状态到领域结论的一对一映射

**问题诊断（Problem Diagnosis）**:
GeoRBF 根据 validation、preflight、certificate 与 canonical 复验证据对问题作出的语义结论，与任何单次求解尝试终止状态分别记录。
_Avoid_: solver enum、日志字符串、未经复验的后端结论

**版本化数值政策（Versioned Numerical Policy）**:
共同标识 scaling、内部表示、rank 语义、分解验收、精度升级、恢复与查询可靠性边界的不可变审计契约；改变其中任一成功或失败判定必须使用新 id。v0.2.1 只执行 `georbf-v2`；`georbf-v1` 保留为 v0.2.0 历史身份，不得静默改写或作为缺陷兼容路径重新启用。
_Avoid_: backend 版本、可变容差集合、同 id 改语义、crate patch 版本即政策 id

**v0.2.1 compatibility boundary**:
保持 v0.2.0 public signatures、领域关系语义、物理验收包络与真实诊断结论；允许 `georbf-v2` 修正旧 relative-spectrum 误判并改变容差内浮点结果，不承诺 bitwise identity。新增证据 accessor 与 non-exhaustive error variant 不要求调用方新增配置。
_Avoid_: 静默 policy 变更、破坏性 API 删除、bitwise promise、保留已知误诊、改写历史 validation

**Canonical 物理验收包络（Canonical Physical Acceptance Envelope）**:
逐关系应用的 `1e-10 * characteristic_scale + 1e-8 * relation_reference_scale` 物理容差公式；`georbf-v2` 保持该领域成功语义不变。内部 basis、factorization、solver 与 query 误差界彼此独立，不能扩大此包络。
_Avoid_: condition-adjusted tolerance、fixture 特例、聚合 residual tolerance、backend feasibility tolerance

**后端契约违反（Backend Contract Violation）**:
后端声称给出可接受候选，但恢复后的后端标准形式不满足该适配器所声明数值契约的失败。
_Avoid_: 恢复验收失败、近似成功、用户输入无效

**问题正则化（Problem Regularization）**:
通过修改 kernel pairing、FieldEnergy、soft objective 或 hard feasible set 来改善数值行为的 ridge、jitter 或 penalty；GeoRBF 禁止隐式问题正则化。
_Avoid_: 后端分解正则化、可逆换基、显式领域 soft loss

**后端分解正则化（Backend Factorization Regularization）**:
后端仅为其内部 KKT 分解采用、由版本化政策完整记录且不写回 Canonical Problem 或 solver-independent form 的数值稳定化；它只可产生候选，不能证明 rank、可行性或成功，候选仍须对未正则形式完整 Recover and Verify。
_Avoid_: kernel ridge、matrix jitter、隐藏设置、正则后 residual 即成功、hard 转 soft

**恢复验收失败（Recovery Verification Failure）**:
候选满足后端标准形式契约，但经 recovery map 返回 Canonical Problem IR 后未通过物理语义、round-trip 或 provenance 验收的失败。
_Avoid_: 后端契约违反、放宽容差后的成功、不可行问题

**直接输入冲突（Direct Input Conflict）**:
在 lowering 前即可由局部代数、几何或关系图证明的 hard 领域关系矛盾。
_Avoid_: 欠约束、未经证明的求解失败、soft 残差

**不可行问题（Infeasible Problem）**:
全部单条输入均合法、但完整 hard affine/conic 可行集经可靠判定为空的问题。
_Avoid_: Rank Deficient、数值失败、迭代未收敛、soft 冲突

**来源可定位冲突见证（Source-Localized Conflict Witness）**:
能够从求解或代数证书完整恢复到原始 canonical hard 关系、列出参与矛盾的 SourceId，并以 residual 与分离裕量复验不可同时满足性的充分见证。它可以确定性删减，但不承诺全局最小基数。
_Avoid_: backend infeasible 状态、无来源证书、数值失败、仅列出 presolve 行、全局最小冲突集承诺

**可解释秩亏（Interpretable Rank Deficiency）**:
由结构或可重构代数证据证明存在非零 canonical functional 组合、polynomial mode 或 quotient field mode 实际作用为零，并能恢复该依赖来源的秩缺失。原始或缩放矩阵的条件数、相对小谱值、近零 pivot 或 Cholesky 失败本身都不是秩亏证明。
_Avoid_: 病态矩阵、小正谱值、数值判定灰区、不可行问题

**正有效场模态（Positive Effective Field Mode）**:
在 Cubic 的完整 \(\Pi_1\) quotient 中产生非零场变化且 FieldEnergy 严格为正的模态；一旦可靠地区分于零，它相对最大模态可以任意小，仍必须保留。
_Avoid_: 数值噪声、小谱值截断、一致冗余、ridge 保留的伪模态

**数值判定灰区（Numerical Decision Gray Zone）**:
主精度证据及版本化政策要求的有界精度升级仍无法可靠证明相关模态为正、代数为零或具有负曲率的问题状态；它产生数值失败而不是确定的秩亏或不可行诊断。
_Avoid_: 可解释秩亏、不可行问题、后端声称成功

**派生 KKT 风险证据（Derived KKT Risk Evidence）**:
solver-facing KKT 的 condition estimate、pivot、inertia、refinement 与 backward error；它们指导 attempt 与精度升级，但不能脱离 canonical algebraic recovery 单独证明 Rank Deficient、Infeasible 或成功。
_Avoid_: KKT condition gate、relative-spectrum rank、backend inertia 即领域诊断、忽略 exact zero assembly row

**拟合问题规模**:
领域观测完成 lowering、但尚未经过求解器 presolve 时形成的标量约束维数；中心系数、辅助变量、锥块和实际 KKT 维数是同时记录的独立规模指标。
_Avoid_: 观测数量、采样点数量、presolve 后约束数

**Dense Hermite 500 release boundary**:
由未修改的 `same_horizon_500.csv` 定义的 v0.2.1 强制容量边界：500 条 fixture records、1,000 个观测 SourceId、一个 gauge SourceId、2,000 个 field representers、1,996 个 Cubic quotient 模态和 2,501 条 canonical hard scalar relations；当前未压缩求解形式约有 2,001 条 hard rows。允许的可重构压缩可以改变 solver row 数，不能改变其余 canonical 计数或全来源参与。
_Avoid_: 500-constraint case、抽样性能集、2,000-constraint 泛称、仅验证 fixture

**Dense Hermite 500 来源映射**:
每条 fixture record 产生一个 Horizon member SourceId 和一个 Directed Normal SourceId，并以共同的稳定 row 前缀关联；独立 gauge 使用第三类 SourceId。参与账本按 500 条 records 汇总，关系恢复与冲突诊断保留精确 channel SourceId。
_Avoid_: 复用 SourceId、只给 normal 编号、虚构复合公开观测、把 gauge 归给某一数据行

**Dense Hermite 500 problem profile**:
500 个位置属于一个 hard Horizon；每个单位法向产生 hard Normal Direction 与显式 hard `MinimumNormalSlope = 1.0` field-unit/input-length；该 Horizon 的 additive gauge 为 `0.0`。单位法向只表达方向，不是幅值已知的完整 Gradient。
_Avoid_: GradientObservation、epsilon slope、soft normal、无 gauge、从数据估计 slope

**Dense Hermite 500 数据身份**:
固定回归 fixture 的规范 LF 字节流、行序和数值文本；其 SHA-256 为 `dec2c30c361bc5341d42838e5d72e13e9372058997e288fa81ca482246e1c1a0`，包含一个表头与 500 条有限位置及单位法向记录。换行转换、重新排序、舍入或数值等价改写都不是同一回归身份。
_Avoid_: CRLF 摘要、解析后重写、容差内修改、生成式 fixture、抽样副本

**Dense Hermite 500 validation profile**:
所有平台校验 fixture 身份与 canonical 计数；标准 profile 在 Linux x86_64 优化构建中执行完整拟合，release/tag profile 在五个受支持原生目标上执行同一完整拟合。平台间比较 canonical 成功语义与物理验收包络，不比较系数或查询结果的 bitwise identity。
_Avoid_: debug-only dense fit、单平台 release、平台特定抽样、bitwise 浮点门槛

**Dense Hermite 500 acceptance bundle**:
要求 `georbf-v2` 成功保留 2,000 个 representers、rank-4 完整 \(\Pi_1\) 与 1,996 个 quotient 模态，不截断、不做问题正则化，并对 2,501 条 canonical hard scalar relations 及最终模型在全部 500 个位置的 value/gradient/normal/slope 查询执行全来源物理验收。fixture 身份、FieldEnergy、objective、side condition、basis/query round-trip 和 Representation Evidence 任一失败都阻止发布。
_Avoid_: fit status only、部分 SourceId、只查训练矩阵、跳过最终模型查询、fixture 特例容差

**Dense Cubic Hermite representation pathology**:
v0.2.0 将 canonical representer coordinates 的病态程度误作场空间代数秩，并在无正则 Cholesky 已成功后仍以全局相对谱阈值拒绝严格为正的局部分辨率模态。显式 nullspace 逐列物化的额外 \(O(n^3)\)、representation failure 证据丢失和 quotient dimension 猜测错误是伴随缺陷。
_Avoid_: 坏数据、重复观测、hard 不可行、Clarabel 失败、Surfe 差异

**结构性性能边界（Structural Performance Boundary）**:
在尚无 manifest-pinned 专用 runner 时，以算法复杂度、checked peak capacity、正常路径分解次数和 release workload 完成性定义的性能契约；普通 CI 的 wall time 与资源趋势只作回归证据，不构成跨机器 SLA。
_Avoid_: GitHub runner 秒数 SLA、只测小样本、成功后全谱 SVD、未规划内存、查询二维矩阵

**逻辑查询批量**:
作为一次求值操作提交的有序查询位置集合；内部自动分块或流式处理不改变其输入顺序及逐点等价语义。
_Avoid_: 查询—中心矩阵、内部分块

**问题快照**:
一次批量拟合所依据的完整原始输入与审计记录，包括观测、约束、核参数、各向异性及求解配置；任何影响解的变化都会形成新的快照。它保留调用方陈述，不承担规范化后的计算语义。
_Avoid_: 可变问题、增量状态、规范计算表示

**Canonical 参与（Canonical Participation）**:
一个原始来源的完整 hard 关系属于规范可行集，或其完整 soft 残差与损失属于规范目标；它不要求来源与求解行或场表示基一一对应。
_Avoid_: 求解行参与、representer 参与、仅验证观测、抽样参与

**全来源恢复验证（All-Source Recover and Verify）**:
把候选解恢复到物理场后，对每条原始 canonical hard 关系逐项验收，并恢复每条 soft 关系的残差与损失贡献；所有 SourceId 都必须出现在参与和验证账本中。聚合 residual、压缩行验证或 representer 子集验证不能替代它。
_Avoid_: KKT residual 即成功、抽样验证、仅验证保留行、条件数放宽容差

**内部场表示基（Internal Field Representation Basis）**:
与规范关系分离、覆盖同一可容许场空间的内部数学坐标；其数量和顺序可以不同，但转换必须可重构并保留全部规范 hard 与 soft 语义。
_Avoid_: Canonical Constraint、观测子集、采样中心、小模态截断

**Cubic 商空间表示层（Cubic Quotient Representation Layer）**:
位于 Canonical Problem IR 与求解形式之间的派生表示层；它显式保留完整 \(\Pi_1\) 多项式空间，在满足多项式旁条件的 Cubic quotient span 中构造内部场表示基，并维护 canonical representer、内部坐标和逐来源响应之间的可重构映射。
_Avoid_: rank threshold 补丁、原始 representer 矩阵即规范语义、多项式消去后遗失、不可恢复的低秩近似

**Cubic Quotient Representation Module**:
位于 Canonical IR 与 solver-form assembly 之间的纯进程内深 Module；其 Interface 仅负责从完整 fitting functionals `build` 表示、为 canonical functional 生成 `response`、以及把 field coordinates `recover` 为已验证查询表示和 Representation Evidence。Householder、LLT、double-double、capacity 与换基证据全部属于其 Implementation，不建立单 Adapter 的假 Seam。
_Avoid_: shallow factor wrappers、backend-owned basis、公开矩阵 Interface、canonical hard/soft 语义下沉、每个算法一个 Module

**Canonical Solver-Form Assembly Module**:
在 canonical hard/soft 语义与 KKT/QP Adapter 之间维护关系响应、可行集或目标恒等压缩、SourceId recovery graph 和冲突证书的 Module；它只通过 Cubic Quotient Representation Module 的 `response` Interface 获取场坐标，不拥有 field basis。KKT 与 QP 必须消费同一 solver-independent form。
_Avoid_: representation-owned hardness、backend-owned provenance、hard/soft 共用压缩规则、presolve row 即 canonical relation

**Representation Evidence bundle**:
成功与失败共享的 Cubic 表示审计记录，直接报告 canonical、representer、\(\Pi_1\)、quotient 与 recovery 的实际维数，以及换基、无正则分解、精度升级、重构缺陷、来源覆盖和失败阶段。维数不得由其他 row count 推算，condition/spectrum 只作为风险摘要。
_Avoid_: success-only evidence、空 representation failure、猜测 quotient dimension、完整谱即 rank 结论、backend log

**\(\Pi_1\) 可识别性（\(\Pi_1\) Unisolvency）**:
全部 representers 对完整 \(\mathrm{span}\{1,x,y,z\}\) 的 pairing 具有数学 rank 4；小而非零的 affine mode 必须通过稳定换基或精度升级保留。只有证明非零 affine polynomial 被全部 representers 湮灭时才诊断 polynomial rank deficiency，并恢复其物理系数与来源覆盖。
_Avoid_: 自动降为 \(\Pi_0\)、删除线性项、relative-SVD 截断、用 gauge 掩盖模式

**隐式商空间合同变换（Implicit Quotient Congruence）**:
以多项式 pairing 的正交补隐式表示完整 quotient，并用合同变换把 kernel pairing 限制到该空间；正交补不作为观测选择器，trailing quotient block 的每一维都保留。v0.2.1 的首选构造使用 Householder 正交变换。
_Avoid_: 抽取 representer、显式稠密 nullspace 逐列物化、低秩近似、Schur 删除观测

**能量正交商空间基（Energy-Orthonormal Quotient Basis）**:
覆盖完整 Cubic quotient span 的满秩内部场表示基；它通过未经正则化且经过验证的可逆换基，使 solver-facing 场坐标的 FieldEnergy 为欧氏平方范数，同时保留全部正有效场模态。
_Avoid_: FieldEnergy Normalization、谱截断、ridge whitening、jitter、近似低秩基

**经验证的无正则分解（Verified Unregularized Factorization）**:
在不添加 ridge、jitter 或截断模态的前提下保持预期完整维数，并以正模态证据、尺度感知的 backward residual 及可重构换基共同验收的分解。对 full quotient Gram matrix \(G\) 与计算所得 \(\hat L\)，v2 的 LLT 证书必须满足
\[
\eta_G=
\frac{\lVert G-\hat L\hat L^T\rVert_\infty}
{\lVert G\rVert_\infty+\lVert |\hat L||\hat L|^T\rVert_\infty}
\le 10^{-11},
\]
且每个 pivot 的 outward-rounded 区间下界严格大于零；Householder 正交性、\(\Pi_1\) side condition、能量恒等式与 canonical response 往返误差也都必须不超过 \(10^{-11}\)。任一 f64 pivot 区间跨零时只对相应歧义模态启动 Targeted double-double rescue；升级后仍跨零则返回 Numerical Decision Gray Zone。只有证明代数零模态时才可返回 Rank Deficient，不得以 \(\sigma_{\min}/\sigma_{\max}\) 或其他 condition gate 替代该证书。例程返回成功、结果有限或条件估计良好都不能单独构成成功证据。
_Avoid_: LU 成功、Cholesky 成功、有限权重、condition gate、隐式 stabilization

**精度升级（Precision Escalation）**:
当主算术精度无法可靠确认 quotient 模态或分解证据时，在同一 Canonical Problem 上以更高精度重算歧义子空间并重新验证；它不改变目标、约束、维数或容差。v0.2.1 在返回 Numerical Decision Gray Zone 前必须执行一次有界升级。
_Avoid_: ridge、jitter、放宽容差、hard 转 soft、全问题无界高精度重算

**Targeted double-double rescue**:
通过 deterministic symmetric pivoting 隔离 quotient trailing Schur 歧义子空间，并以约 106-bit 的纯 Rust double-double 算术从原始 canonical Cubic pairings 重算该子空间的精度升级；确认的正模态重新接入完整分解，最终 f64 换基仍须全量验证。`georbf-v2` 每次 fit 最多升级 `min(64, quotient_dimension)` 个歧义模态，超限返回灰区而不截断。
_Avoid_: 对已舍入 \(G\) 重复分解、完整稠密高精度 solve、截断 trailing block、native multiprecision 依赖

**Precision-rescue oracle corpus**:
以独立至少 120-digit 参考值验证 double-double arithmetic、Cubic jet/pairing、Schur accumulation，以及小正、真零、负曲率和 64/65 模态边界的分层语料；它必须证明 rescue 保留正模态且不会把零或负模态修成可接受。
_Avoid_: self-oracle、只测最终 fit、只有正例、已舍入矩阵参考、跳过上限边界

**硬约束（Hard Constraint）**:
属于拟合可行集、并要求成功模型在数值验收容差内满足的单条领域关系。
_Avoid_: 大权重、无限权重、自动放松、全模型 hard mode

**软关系（Soft Relation）**:
允许违反、并把显式违反量通过指定损失纳入拟合目标的单条领域关系；其残差保留原关系单位。
_Avoid_: 放宽数值容差、隐式默认权重、全模型 soft mode

**目标恒等 soft 压缩（Objective-Identical Soft Compression）**:
对每个候选场都严格保持全部原始 soft 残差、multiplicity、权重与损失之和相同，并能恢复逐 SourceId 残差及损失贡献的代数压缩。仅有 soft 行空间相同不足以准入。
_Avoid_: 删除重复 soft 观测、合并后改变权重、仅验证 soft 点、二次近似非二次损失

**软边界（Soft Bound）**:
通过上下侧各自独立的非负违反量与显式惩罚软化单侧或双侧仿射界限。
_Avoid_: 共用双侧 slack、标准差、删失或区间似然

**一致冗余 hard 等式（Consistently Redundant Hard Equality）**:
其完整 canonical 仿射行可由保留 hard 等式通过同一恢复映射同时重构左侧算子与右侧目标，因而代数压缩不改变可行集的关系；完全相同的 hard 行只是其特例。压缩后仍保留每个原始 SourceId、恢复系数及逐关系验证。
_Avoid_: 仅左侧相关、近似相关、soft 重复证据、丢弃来源、坐标扰动

**近邻独立观测（Near-Support Independent Observation）**:
support 坐标不同、因而仍定义独立 canonical functional 与场分辨率模态的观测；间距小或法向相同不构成冗余。它只能经一般可逆换基参与稳定表示，不能由 merge radius、snapping 或 tolerance dedup 合并。
_Avoid_: 几何重复、噪声点、可删除小模态、坐标聚类

**代数语义回归语料（Algebraic Semantics Regression Corpus）**:
独立验证一致 hard 冗余、冲突 hard 依赖、近邻独立 supports、相同法向不同 supports、soft multiplicity 及 hard/soft 同左侧共存的固定微型案例；它们证明压缩和诊断规则，但不能替代 Dense Hermite 500 完整回归。
_Avoid_: 主 fixture 抽样、仅成功案例、只检查 solver row count、soft 去重

**非规范对照（Non-Normative Comparison）**:
只用于理解其他实现行为、但不参与 GeoRBF pass/fail、rank、可行性或成功结论的外部结果。Surfe 属于此类；其 LU 返回或有限权重既不是成功证据，失败也不是不可行证据。
_Avoid_: oracle、backend dependency、兼容目标、宽松成功判据

**统计不确定性（Statistical Uncertainty）**:
以标准差或协方差描述观测误差的统计语义，与仅表达优化取舍的惩罚权重不同。
_Avoid_: 任意权重、置信度标签、求解容差

**协方差组（Covariance Group）**:
以有限、对称、严格正定的协方差共同描述一组同量纲软残差的命名统计单元；其交叉项只产生规范的组级总目标贡献。
_Avoid_: 奇异协方差、hard/soft 混合组、虚构的唯一逐观测贡献

**平方惩罚（Quadratic Penalty）**:
以有限正权重惩罚标量或同量纲欧氏向量的平方残差，权重不具有统计不确定性含义。
_Avoid_: 默认权重、零权重、逐坐标 L1、方差

**线性违反惩罚（Linear Violation Penalty）**:
以有限正权重惩罚显式非负标量违反量的线性损失，权重不具有统计不确定性含义。
_Avoid_: 向量分量 L1、默认权重、零权重、置信度

**已求解模型**:
批量拟合成功后产生并永久对应其问题快照的不可变模型，可安全地重复及并发只读查询。
_Avoid_: 在线模型、原地更新模型

**已验证查询表示（Verified Query Representation）**:
从稳定求解表示一次性恢复、并经全来源响应与 FieldEnergy round-trip 验收的 generalized-RBF 与完整 \(\Pi_1\) 展开；它表示同一个已求解场，使查询可以流式扫描 representers，而无需逐查询重做 quotient triangular solve。
_Avoid_: 第二次拟合、查询时换模型、未验证原始系数、每点 \(O(n^2)\) 恢复

**自适应可靠查询（Adaptive Reliable Query）**:
对已验证查询表示以确定性补偿累加计算 value 与 gradient，并依据 `1e-12 * field_scale + 1e-11 * sample_reference_scale` 的分量级物理误差包络只在必要时升级累加精度的 \(O(n)\) 查询；升级始终评估同一个已求解场。升级后仍不能可靠判定时返回结构化 numerical-indeterminate 查询错误。
_Avoid_: 逐项朴素累加、切换模型、查询时重新拟合、每点 quotient solve、静默返回不可靠有限值

**Kernel 契约（Kernel Contract）**:
对一个准入 kernel 不可拆分地规定精确公式、空间尺度、正定性或条件正定性、多项式要求、可用 jet、各向异性协变、FieldEnergy 及准入证据的数学契约；名称本身不构成 kernel。
_Avoid_: kernel 名称或枚举、solver 私有配置、只有求值函数的自定义 kernel

**Kernel 空间尺度（Kernel Spatial Scale）**:
由调用方以输入坐标框架的长度单位显式提供、用于形成无量纲径向距离的严格正长度；它不得从数据或内部坐标归一化静默估计。
_Avoid_: shape、epsilon、自动点间距、FieldEnergy normalization、Global Anisotropy Metric

**支撑半径（Support Radius）**:
紧支撑 kernel 在各向异性度量距离中从非零分支精确进入零分支的严格正物理长度。
_Avoid_: 普通 length scale、近似截断阈值、稀疏求解容差

**场能量（FieldEnergy）**:
由 Kernel Contract 唯一定义的场复杂度度量；SPD kernel 使用 native-space norm，Cubic 使用完整 \(\Pi_1\) quotient 上的 seminorm。
_Avoid_: ridge、jitter、solver 正则化、soft penalty

**场能量归一化（FieldEnergy Normalization）**:
把 kernel 的原生场能量转换为无量纲目标项的显式有限正系数；它决定场能量与 soft losses 的相对强度，并随物理单位协变。
_Avoid_: 隐藏正则化强度、kernel 长度尺度、自动 normalization、solver scaling
