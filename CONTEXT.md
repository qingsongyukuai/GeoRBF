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

**直接输入冲突（Direct Input Conflict）**:
在 lowering 前即可由局部代数、几何或关系图证明的 hard 领域关系矛盾。
_Avoid_: 欠约束、未经证明的求解失败、soft 残差

**不可行问题（Infeasible Problem）**:
全部单条输入均合法、但完整 hard affine/conic 可行集经可靠判定为空的问题。
_Avoid_: Rank Deficient、数值失败、迭代未收敛、soft 冲突

**拟合问题规模**:
领域观测完成 lowering、但尚未经过求解器 presolve 时形成的标量约束维数；中心系数、辅助变量、锥块和实际 KKT 维数是同时记录的独立规模指标。
_Avoid_: 观测数量、采样点数量、presolve 后约束数

**逻辑查询批量**:
作为一次求值操作提交的有序查询位置集合；内部自动分块或流式处理不改变其输入顺序及逐点等价语义。
_Avoid_: 查询—中心矩阵、内部分块

**问题快照**:
一次批量拟合所依据的完整原始输入与审计记录，包括观测、约束、核参数、各向异性及求解配置；任何影响解的变化都会形成新的快照。它保留调用方陈述，不承担规范化后的计算语义。
_Avoid_: 可变问题、增量状态、规范计算表示

**硬约束（Hard Constraint）**:
属于拟合可行集、并要求成功模型在数值验收容差内满足的单条领域关系。
_Avoid_: 大权重、无限权重、自动放松、全模型 hard mode

**软关系（Soft Relation）**:
允许违反、并把显式违反量通过指定损失纳入拟合目标的单条领域关系；其残差保留原关系单位。
_Avoid_: 放宽数值容差、隐式默认权重、全模型 soft mode

**软边界（Soft Bound）**:
通过上下侧各自独立的非负违反量与显式惩罚软化单侧或双侧仿射界限。
_Avoid_: 共用双侧 slack、标准差、删失或区间似然

**冗余关系（Redundant Relation）**:
规范化后与已有 hard 关系完全相同、可为数值求解合并但仍保留全部来源身份的关系。
_Avoid_: soft 重复证据、丢弃来源、坐标扰动

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
