# Surfe 兼容与错误策略

## 兼容目标

GeoRBF 的目标是对冻结 Surfe 非可视化核心的实际可达行为等价，而不是复刻 C++ ABI、
对象所有权、控制台输出或内存模型。来源身份固定为
`290dbe0ab344f4258a4935f05cad0f153f0f69a4`；所有来源引用都必须使用
`path@commit` 形式。

兼容证据按以下顺序组合，而不是任选一项：

1. 冻结源码的定义、分支、布局和调用点；
2. 固定提交 oracle 的可观察结果和中间证据；
3. 解析恒等式、符号/对称性检查和有限差分；
4. `PLAN.md` 中任务专属完成门槛。

文档说明只能帮助解释源码，不能覆盖相反的源码或 oracle 证据。T03 之前没有 oracle，
T04 之前没有正式 fixture/tolerance；因此 T00 只冻结策略，不声明数值 parity。

## 有效输入与离散行为

“有效输入”是同时满足冻结公开入口的形状/状态要求、相关模型约束前置条件、有限数值
要求，以及不会触发 C++ 未定义行为的数据。后续任务必须用源码与 oracle 为每个入口
收窄这一定义，不能用“Rust 更方便”扩大接受集合。

以下离散行为必须精确比较：

- 参数默认值、枚举和兼容名称；
- 配置的接受/拒绝、调用顺序和成功/失败类别；
- 约束排序、`1e-3` 同点判定、分类去重和精确 level 分组；
- reference point、increment pair、索引、自由度、矩阵分区和项顺序；
- equality/inequality 符号、solver 分支、active/reconstruction 选择；
- 对有效输入可达或不可达的模型/功能。

不得用容差比较离散字段，不得额外接受大小写变体、别名或自动修复的配置。字符串和
默认值的最终表由 T06/T30 固定。

### T06 参数、名称与异常映射

T06 对 `surfe_lib/modelling_parameters.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
和 `surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 做了以下离散固定：

- `Parameter_Types` 的七组枚举保持声明顺序；`Parameters`、`InternalParameters` 和
  `InputParameters` 保持逐字段默认值。Rust 使用安全的值类型和拥有所有权的字符串，
  不复制 C++ ABI 或内存布局。
- `SetRBFKernel(const char*)` 只接受源码中的九个逐字节精确名称：`r3`、`WendlandC2`、
  `r`、`Gaussian`、`Multiquadratics`、`Multiquadratics3`、`Thin Plate Spline`、
  `Inverse Multiquadratics`、`MaternC4`。不修剪空白、不忽略大小写、不增加别名。
- 公开 `Surfe_API(int)` 的模型代码固定为 1 Single Surface、2 Lajaunie、3 Vector
  Field、4 Stratigraphic Horizons、5 Continuous Property；这与 `ModelType` 中后三项的
  声明顺序不同。GeoRBF 的模型文本形式只使用冻结 C++ 枚举标识符并精确往返；冻结
  Surfe 本身没有接受模型字符串的公开入口，因此这不是对 Surfe 输入集合的扩张。
- `SetRegressionSmoothing` 和 `SetGreedyAlgorithm` 忽略传入布尔值并无条件把对应字段写为
  `true` 的源码行为由参数级兼容测试固定；这不表示 Greedy 调用链已可达，后者仍归 T31。

T06 还逐项映射了
`surfe_lib/grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 23 个具体
异常类。Rust 代码以 `Error` 枚举类别匹配；C++ 类名、原始 `what()` 文本和
`SurfeExceptions` 的嵌套文本仅作为可追踪诊断证据，不作为程序化类别。冻结头文件使用
`std::string` 却未自行包含 `<string>`；oracle 探针按冻结构建的包含顺序补齐标准头，未修改
reference，也未把这一编译卫生缺陷复制到 Rust。

### T07 约束几何与安全输入边界

T07 完整核对了
`surfe_lib/modelling_input.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
的 `Point`、`Interface`、`Inequality`、`Planar`、`Tangent` 构造与方向方法，并补充核对
`surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 azimuth 入口：

- 坐标、`c`、level、scalar/vector field 均由拥有所有权且初始化的 Rust 值保存；不复制
  冻结默认构造函数的未初始化字段。Planar normal bounds 和 Tangent angle bounds 在
  C++ 中设置前未初始化，Rust 用 `Option` 明确表示“尚未计算”。
- normal 到方向的换算保持 `dip=acos(nz)`、`dip_direction=atan2(ny,nx)`、负方位角加
  360°、`strike=360°-dip_direction`；只有 `nz < 0` 是 overturned。strike/dip 到 normal、
  dip/strike vector 和四角 normal envelope 逐运算保留冻结公式与度单位。
- azimuth 是 dip direction；`azimuth >= 90°` 时 `strike=azimuth-90°`，否则
  `strike=azimuth+270°`。公开文档限定 polarity 只能为 0 upright 或 1 overturned，Rust
  用离散枚举拒绝其他整数。
- 有效的非单位 normal/tangent 分量按原值保存，不做隐式归一化；只有角度构造产生的
  normal 和 `getDipVector` 依冻结源码归一化。这样不会以“修复”改变传入方向对后续矩阵
  的尺度影响。
- 冻结 C++ 探针证实零 normal、零 tangent 和非有限坐标会被接受，`nz > 1` 会令 dip
  成为 NaN。它们不属于有效输入；Rust 分别返回 `ZeroNormal`、`ZeroTangent`、
  `NonFiniteInput` 或 `NormalZOutOfRange`，不把 NaN/未定义状态带入后续数学路径。

T07 不实现 `Point` 排序/同点判定、空间算法、残差 Greedy 行为、reconstruction setter、
矩阵、求解或 fitted model；这些仍由 T08/T09/T21/T29/T31 的既有任务关闭。

### T08 排序、同点、去重与 level/reference 分组

T08 核对了
`surfe_lib/modelling_input.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、
`surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 和
`math_lib/math_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的实际离散语义：

- `Point::operator<` 只按 x、y、z 做逐坐标精确字典序，不读取 `c` 或约束 payload；
  `-0.0` 与 `+0.0` 比较相等。冻结源码不存在计划提示的 `Point::operator==`，实际同点
  判定是自由函数 `collocated`：三个轴的绝对差必须各自严格小于 `Epilson=1e-3`，恰在
  边界即不是同点。
- 四类约束各自排序后各自执行相邻 `std::unique(collocated)`；绝不跨类别删除，也不把
  非传递的同点关系改写成全局聚类。Rust 使用相同的“与最近保留项比较”语义，保留
  排序后每段的第一项及其 level/normal/tangent payload。
- C++ `std::sort` 没有规定比较等价项之间的顺序。冻结 g++ oracle 对 T08 重复矩阵保留
  了输入中首个精确等价项；Rust 的稳定排序适配将这一结果变成跨平台确定规则，同时不
  改变任何不同 `(x,y,z)` 点的源码顺序。冲突 payload 的精确同坐标重复不再依赖标准库
  私有排序实现。
- interface 和 inequality level 都只按 `==` 精确去重并从大到小排列，不使用容差；相邻
  binary64 level 仍属于不同组。每个 interface level 的 reference 是清洗后位置顺序中的
  首点。singleton level 仍保留 level/reference，但依冻结 `_get_interface_points` 从多点
  组列表移除；同层自由度前置计数为每个保留多点组的 `len-1` 之和。
- `Math_methods::sort_vector_w_index` 不是普通稳定排序。Rust 保留其 7 项 insertion
  threshold、partition 调度、50 项显式栈、重复值和 signed-zero 的 value/index 联动顺序；
  长度不匹配与冻结函数一样返回失败且不改输入。

T08 不建立空间距离/最近邻、increment pair、模型 layout 或矩阵；这些仍分别属于
T09、T15/T25/T26、T16 和 T17。

### T09 空间辅助算法

T09 完整核对了
`surfe_lib/modelling_input.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
的约束转点、距离、最近/最远邻、平均最近邻、`spatial_metrics`、bounds、最远点对、
目标距离点、极值采样、最大距离、轴向变化顺序和四类约束平均距离：

- `distance_btw_pts` 是包含 `c` 坐标的四维欧氏距离；bounds 和极值轴仍只使用 x/y/z。
  约束转点顺序精确为 inequality、interface、planar、tangent，并保留各点的 `c`。
- 单最近邻跳过所有距离精确等于零的候选，以严格 `<` 保留首个并列索引；多最近邻先
  过滤零距离，再复用 T08 的冻结 indexed sort 顺序。冻结函数在请求数仍包含已过滤点时
  可越界；Rust 只返回实际存在的最多 `n` 个索引。空候选或没有非零邻居分别返回
  `EmptyPointSet`/`NoNonzeroNeighbour`，不传播 `-1` 哨兵。
- 最远邻、跨集合最远索引、目标距离最近索引和最远点对都保留严格比较产生的首个并列
  索引。两个精确重复点的最远点对按冻结双循环仍是 `[0,0]`；少于两点、空集合或非有限
  目标距离使用类型化 `SpatialError`，不复制 `pts[0]` 越界或无效索引。
- `avg_nn_distance` 按“不同索引”而不是“非零距离”选最近点，因此重复点的最近距离可为
  零；空集和单点都精确返回零。四类 `Constraints::compute_*_avg_nn_distance` 保持这一
  数值定义。冻结 `compute_avg_nn_distances` 把四值写入可变缓存；Rust 返回拥有所有权的
  `ConstraintAverageNearestNeighbourDistances` 快照，消除隐藏共享突变而不改变计算结果。
- `spatial_metrics` 先按 T08 规则排序和相邻去重，再以四维距离求平均最近邻并除以二，
  同时输出 x/y/z bounds。冻结空输入仍返回 `true`、零 resolution 以及
  `DBL_MAX/-DBL_MAX` bounds；这些哨兵不构成可用空间参数，Rust 对空输入返回
  `EmptyPointSet`。单点仍成功，resolution 为零且每轴上下界相等。
- extremal 索引和轴向变化顺序继续使用冻结 `sort_vector_w_index` 的并列顺序；各轴 range
  从大到小，range 完全相等时顺序为 Z、Y、X。空 extremal 输入在 C++ 中会访问越界，
  Rust 安全拒绝；最大两点距离则与冻结源码一样为空集/单点返回零。

同一文件的四个 `Get_*_STL_Vector_Indices_With_Large_Residuals` 消费本任务的距离和平均
距离，但属于 Greedy 的候选选择、残差阈值与追加策略；依固定计划仍归 T31，本任务没有
提前移植。`continuous_property.cpp` 中注释掉的 extremal 调用也未被误报为可达能力。

## T11 Lagrangian/unisolvent 选择和退化处理

T11 完整核对了
`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
的 `Lagrangian_Polynomial_Basis::{_get_unisolvent_subset,_initialize_basis,poly,poly_dx,
poly_dy,poly_dz}` 和
`surfe_lib/grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
的 `failurecreatinglagrangianpolynomialbasis`：

- 选择使用第一个拥有严格最大点数的 interface group；x/y/z 分别复用冻结 indexed sort，
  再从最大 range 轴取两端、从次大 range 轴取排除两端后的两端。GeoRBF 暴露所选 group
  与四个原始组内索引作为离散证据。
- 冻结三个最大轴判断是独立 `if`，不是 `else if`。最大 range 精确并列时会累计多于四点
  并抛出 Lagrangian 构造错误；即使四个输入点在数学上本可 unisolvent，Rust 仍保留该
  有效可观察的拒绝分支。
- 若初选四点全落在 x、y 或 z 常量平面，冻结路径先按原输入顺序尝试替换最后一点；仍
  找不到不同坐标时，对所选副本的第一个点在对应轴加精确 `Epilson = 1e-3`。Rust 保留
  这一合成点、原始索引和后续系数结果，不修改调用方的约束。
- 一阶四函数的 `[constant,x,y,z]` 系数、求值运算顺序和常量导数按源码显式行列式公式
  移植；level 与 `c` 不进入 polynomial，且任何非零、有限行列式都会尝试，不按大小或
  条件数预拒绝。`1e-12` 尺度的近共面用例与冻结 probe 仍逐位一致并成功。
- 冻结路径不检查一般斜共面、共线或重复选择的零行列式；probe 实测构造返回后 basis
  含非有限值。此类输入不能形成可用 Lagrangian basis，GeoRBF 不传播 NaN/Inf，而统一
  返回 `Error::LagrangianBasisCreationFailure`。空组、最大组少于四点和并列最大轴则与
  冻结源码本来就抛出的同一错误类别精确一致。

这一安全拒绝只关闭冻结的非有限无效状态；所选点有效且行列式非零时，选择、系数、
Kronecker 恒等式、导数和任意点求值均保持冻结数值语义。Modified Kernel 的消费仍归
T14，本任务没有提前实现核消去组合。

## 数值行为

核、两点导数、混合 Hessian、anisotropy、Modified Kernel、矩阵/RHS、标量场、
梯度场、QP 目标与迭代证据使用 T04 定义的分层容差。T04 之前不得写入临时宽松阈值
作为事实。

权重向量在病态或非唯一系统中不要求逐位相同；必须比较有限性、后验残差、约束
可行性、目标值（适用时）和预测场。不得只因条件数差而预先拒绝系统：先按冻结
Surfe 的分支尝试求解，再以有限性、残差和可行性分类结果。全局 parity 完成前不得
用低秩、稀疏、局部近似、FMM 或不同优化算法改变求解定义。

核及其导数必须同时通过冻结 C++ golden、解析恒等式/符号/对称性和有限差分三角
验证。三者冲突时先保存可复现实证，再记录兼容决定；不得放宽容差或删掉用例。

## 错误分类

Rust 使用稳定、可匹配的类型化错误，不以 C++ 异常消息文本作为类别。T06 将冻结
`grbf_exceptions.h` 的逐项映射；最低层级必须能区分：

- 未拟合、拟合后参数/约束已变化；
- 未知模型、未知核、非法参数或调用顺序；
- 数组/矩阵维数错误、空输入、非有限输入；
- 约束不足、几何退化、Lagrangian/unisolvent 构造失败；
- 矩阵/RHS 验证失败、LU 数值失败；
- QP 不可行、不收敛或数值失败；
- reconstruction 阶段失败；
- 内部不变量破坏。

如果冻结源码对无效输入发生越界、未初始化读取、悬垂访问、数据竞争或其他未定义
行为，GeoRBF 必须安全拒绝并使用最接近的稳定类别。该安全修复必须记录，但不能改变
有效输入的数学结果。

## 源码缺陷处理

发现疑似缺陷时采用固定顺序：

1. 记录准确的 `path@commit`、符号和调用证据；
2. 在 oracle 可用后增加兼容测试，保存有效输入上的可观察行为；
3. 区分有效输入结果、不可达代码和未定义行为；
4. 在本文记录最终处理，再实施安全 Rust 行为；
5. 不因“数学上更合理”静默改变标准有效输入结果。

T00 阅读 `surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
时登记了以下后续核验项；它们目前是源码观察，不是已经证实的兼容结论：

- `Surfe_API::SetTangentConstraints` 调用 `AddPlanarConstraintwNormal`，需要 T02 做
  缺陷/可达性分类，并由 T03/T30 的公共入口证据决定有效输入行为。
- `Surfe_API::Surfe_API(int)` 未显式初始化全部状态布尔值，需要 T02 分类；Rust 不会
  复制未初始化读取，T29/T30 应验证安全状态机。
- `SetRegressionSmoothing` 和 `SetGreedyAlgorithm` 没有把传入的启用布尔值写入参数，
  而是无条件写 `true`；分别归入 T02/T30 和 T02/T31 的实际调用链核验。

这些发现不创建新任务，也不允许 T00 提前实现修复。

### T02 源码分类结果

T02 已在 `symbol-classification.md` 对冻结声明、定义和活动调用集做双向核对。
以下仍是源级事实，必须由 T03 oracle 和各归属任务的兼容用例固定：

- `surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的
  `Surfe_API::SetTangentConstraints`、`Surfe_API(int)`、两个 bool setter 和批量评估路径
  分别显示错误约束目标、未初始化状态读取、忽略入参和 stale-state/数据竞争。
- `surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的
  `GRBF_Modelling_Methods::run_greedy_algorithm` 以及三个模型的
  `convert_modified_kernel_to_rbf_kernel` 都没有冻结公开根调用边；有完整函数体不等于可用能力。
- `surfe_lib/basis.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 中 `R/AR` 的 30 个导数
  符号直接或间接进入整数 `-666` 哨兵；其他 13 个具体核的同名方法定义完整。
- `surfe_lib/modelling_input.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的
  `spatial_metrics` 对空输入仍成功，而 Greedy tangent residual selector 在特定分支漏返回；
  后者的唯一上游不可达。
- `surfe_lib/continuous_property.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 RHS 循环、
  `surfe_lib/stratigraphic_surfaces.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的重复 pair
  生成，以及四模型共享可变 kernel 的向量评估均已标为缺陷；Rust 路径必须
  在不改变有效串行数值定义的前提下消除越界和数据竞争。

T02 没有运行 C++ oracle、生成 golden 或定义新语义；上述项仍分别由
T12/T21/T26/T27/T29/T30/T31 关闭，不创建新任务。

## 状态、批量与并发语义

`Surfe_API::ComputeInterpolant` 的可观察流水线为约束清洗、输入处理、方法参数、basis
设置、solver 设置，再把 interpolant 标为可用。标量和向量单点/批量结果属于兼容
范围；控制台进度文字、刷新频率、OpenMP 调度和不安全计数器不属于范围。

最终 `FittedModel` 应不可变并支持安全并发只读。批量和逐点评估必须在分层容差内
一致且保持输入顺序。线程安全化不得通过改变核、矩阵、求解或预测定义实现。

## 能力声明规则

- 只有定义、调用链和有效 oracle 用例共同证明的功能才标为可达。
- TODO、空实现、仅声明、仅 GUI 可达或仅调试可视化的路径不得宣称支持。
- 正式 parity 必须等 T32 全部通过；性能必须等 T33 同机实测；发布结论必须等 T34。
- reference 缺失、命令未运行、用例跳过或结果不确定时，状态只能是未验证/未通过。
