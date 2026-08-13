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

## T12 isotropic 核、导数和源码冲突

T12 完整核对了
`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
的 `RBFKernel::{radius,basis_pt_pt,basis_pt_planar_x,basis_planar_x_pt,
basis_pt_planar_y,basis_planar_y_pt,basis_pt_planar_z,basis_planar_z_pt,
basis_planar_planar}`，以及九个普通核 `Cubic/Gaussian/MQ/MQ3/TPS/IMQ/R/
WendlandC2/MaternC4::{basis,dx_p1,dx_p2,dy_p1,dy_p2,dz_p1,dz_p2,dxx,dxy,
dxz,dyx,dyy,dyz,dzx,dzy,dzz}`：

- 半径严格包含 `x/y/z/c` 四个差值；公开一阶导数和混合 Hessian 仍只作用于 x/y/z。
  Gaussian 对参数平方，MQ/MQ3/IMQ 将参数直接加到 `r²`，WendlandC2 将参数作为
  cutoff，MaternC4 将其作为尺度；GeoRBF 不在核层偷偷重解释参数。
- Cubic、Gaussian、MQ、MQ3、TPS、IMQ、WendlandC2 和 MaternC4 的两点一阶导数与
  3×3 mixed Hessian 均按冻结表达式和符号移植。WendlandC2 只在 `radius > cutoff`
  时提前返回正零；恰在支撑边界仍执行源码公式，因此保留其可观察 signed-zero 结果。
- `R::basis` 可达并返回四维半径，但 `R` 的 12 个直接导数体抛整数 `-666`，另外三个
  Hessian 别名转入同一哨兵；冻结源码没有可供 parity 的线性核导数。GeoRBF 保留
  value 能力，并以稳定的 `KernelError::LinearDerivativeUnavailable` 安全替代未被
  `std::exception` 捕获的整数异常，不臆造数学导数或宣称九核导数均已实现。
- TPS 与 WendlandC2 的对角 Hessian 源表达式在 `c1 != c2` 时没有把 `c_delta²` 加入
  所有链式项，因此与对完整四维径向函数做有限差分冲突。公开 Surfe 约束构造的
  `c` 固定为零；GeoRBF 对标准 `c1=c2` 路径同时满足源码、解析恒等式和有限差分，
  对显式不同 `c` 则保留冻结 hexfloat 结果并用兼容测试记录，不静默“修正”源码。

冻结 probe 覆盖九核 separated/zero/near、非零 `c_delta`、Wendland 支撑内/边界/外和
oblique interior；其全部有效结果与 Rust 测试逐位一致。各向异性、Modified Kernel 和
统一 Value/Planar/Tangent 泛函仍分别归 T13、T14、T15，没有在 T12 提前实现。

## T13 全局各向异性和 binary32 中间语义

T13 完整核对了
`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
的 `RBFKernel::{get_global_anisotropy,scaled_radius,_Global_Plunge,_Transform}`，以及
`ACubic/AGaussian/AMQ/ATPS/AIMQ/AR::{basis,dx_p1,dx_p2,dy_p1,dy_p2,dz_p1,
dz_p2,dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz}`；同时核对了
`surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
的 `GRBF_Modelling_Methods::create_rbf_kernel` 各向异性工厂：

- normal 外积的六个和先按 C++ `double` 输入顺序累加，再赋给 `Matrix3f`；GeoRBF 的纯
  Rust 3×3 self-adjoint 路径保留冻结 Eigen 的 binary32 缩放、三对角化、隐式 QR、升序
  特征值/向量配对和每步舍入，不以 binary64 特征分解替代。
- `_Global_Plunge` 是最小特征值对应的第一列特征向量。冻结源码写入后没有任何读取点；
  GeoRBF 只把它作为 parity/debug 证据暴露，不赋予新的模型含义。
- 只有 `eVals(0)` 与 `eVals(1)` 在严格 `< 0.0001f` 时提升到 `0.0001f`；第三项不截断。
  支持矩阵保持 `V * diag(1,sqrt(e1/e0),sqrt(e2/e0)) * V^T` 的 binary32 运算顺序。
  数学上它是对称正半定，但冻结浮点乘法可让镜像项相差一个 f32 ULP；离散布局精确，
  对称/正交性质按 T04 `anisotropy_f32` 层检查。
- `scaled_radius` 只变换 x/y/z 差，不读取 `Point::c`；这与 T12 普通核半径包含 `c` 的
  行为有意不同。后续一阶导数使用 `T^T(T delta)`，混合 Hessian 使用 binary32
  `T` 列内积扩展到 binary64 后的冻结表达式。
- 冻结工厂只有 Cubic、Gaussian、MQ、TPS、IMQ 和 R 六个各向异性类；启用 anisotropy
  后请求 MQ3、WendlandC2 或 MaternC4 会进入 `unknown_rbf`，GeoRBF 明确返回
  `AnisotropyError::UnsupportedKernel`，不臆造三个类。`AR::basis` 可用，而 15 个直接或
  间接导数符号仍进入整数 `-666`；Rust 与 T12 一样返回
  `KernelError::LinearDerivativeUnavailable`。
- 少于两个 planar 与冻结 `failurecomputingglobalanisotropy` 对应。冻结源码不检查
  eigensolver 状态且无条件返回 `true`；对有限 normal 导致 binary32 covariance/特征路径
  非有限或不收敛的无效状态，GeoRBF 安全返回 `NonFiniteComputation` 或
  `EigenSolverFailure`，不传播 NaN/Inf 支持矩阵。有效有限输入的数学结果不变。

冻结 T13 probe 的 oblique case 对 eigenvalues、transform、plunge 逐 binary32 bit 固定；
六个核的值路径和五个完整导数核的值、两点一阶导数、3×3 mixed Hessian 共 81 个
binary64 输出与 Rust 逐位一致。单位协方差极限、退化/`0.0001f` 截断两侧、normal 统一
缩放、支持矩阵对称性质及两层有限差分同时通过。Modified Kernel 和统一线性泛函仍分别
归 T14/T15，T13 没有提前实现。

## T14 Modified Kernel 全组合与固定对角语义

T14 完整核对了
`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
的 `Modified_Kernel::{Modified_Kernel,basis_pt_pt,basis_pt_planar_x,
basis_planar_x_pt,basis_pt_planar_y,basis_planar_y_pt,basis_pt_planar_z,
basis_planar_z_pt,basis_pt_tangent,basis_tangent_pt,basis_planar_planar,
basis_tangent_tangent,basis_planar_tangent,basis_tangent_planar}` 及复制/析构边界，
并补充核对
`surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
的 `GRBF_Modelling_Methods::setup_basis_functions` 构造异常包装：

- GeoRBF 直接复用 T11 的四项 Lagrangian basis 与 T12/T13 的普通/各向异性核。
  value、分别作用于第一/第二点的三项梯度和 row-by-column 3×3 mixed Hessian 均保留
  冻结源码的 `base - t1 - t2 + t3 + t4` 以及 j/k 累加顺序。Hessian 的 `t4`
  特别保留源码的 `p2(k) * K(Uj,Uk) * p1(j)` 乘法次序；交换因子会使部分结果产生
  1–2 个 binary64 ULP 的可观察差异。
- 冻结 `t3` 是 `sum_j p1(j)*p2(j)`，只有 `t4` 的非对角项实际求
  `K(Uj,Uk)`；也就是说源码把四个对角值固定为 1，而不是读取实际
  `K(Uj,Uj)`。GeoRBF 不将它“修正”为完整 Gram matrix。Gaussian shape=1 的
  `K(Uj,Uj)=1` 路径通过了四个 unisolvent 点的双向消去性质；其他核仍严格保留
  冻结固定对角语义。
- Value/Planar/Tangent 的全笛卡尔组合由同一 modified value/gradient/Hessian
  收缩得到，但每个公开入口保持冻结符号和作用点。Tangent 使用传入原始方向尺度，
  不在核层归一化；统一 `Value/Derivative/Tangent/Difference` 泛函仍归 T15。
- 普通核九个 value 路径、八个完整导数路径和各向异性六个 value/五个完整导数路径
  均可作为底层；`R/AR` modified value 仍可用，任一导数组合继续稳定返回
  `KernelError::LinearDerivativeUnavailable`，对应冻结整数 `-666`，不臆造导数。
- `Modified_Kernel` 直接构造失败的内层链仍是 Lagrangian basis 创建失败，而
  `setup_basis_functions` 再包装为 `failurecreatingmodifiedkernel`。GeoRBF 的公开
  `ModifiedKernel::from_isotropic/from_anisotropic` 对应这一安全外层构造边界并返回
  `Error::ModifiedKernelCreationFailure`；T11 的独立构造 API 仍保留
  `LagrangianBasisCreationFailure`。冻结复制构造浅拷贝 `_aLPB` 且析构不释放它，
  会泄漏内存；Rust 以拥有所有权的安全 `Clone` 替代，不复制泄漏或悬垂风险。

冻结 T14 probe 对 Cubic 普通/各向异性各记录 separated 的 16 项 value/gradient/
Hessian、9 项 Tangent 组合和 zero-distance 的 16 项，共 82 个 binary64 输出，Rust
全部逐位一致；Gaussian 固定对角消去、linear value/`-666` 和无效 Lagrangian 构造链
也由同一确定性 probe 记录。交换对称、消去性质和两层中央有限差分同时通过，未装配
完整系统、未生成 T32 正式 fixture，也未宣称全局 parity。

### T17 矩阵、RHS 与 smoothing 装配

- 稠密矩阵按 T16 的同一行列标签以 row-major binary64 存储；约束块由 T15 泛函逐项
  求值。Difference/Difference 保留冻结 `(v1-v2)-(v3-v4)` 括号，polynomial block
  保留源码的 `P/PT/0`、完整或截去常数项的顺序和 tangent 三项乘加次序。
- 普通 QP 保留 Single Surface inequality level 的逐行正负号、Stratigraphic 的前缀
  inequality/后缀 equality 切片和各自 RHS。restricted-range 保留
  `b <= A*x <= b+r`，并在装配前以冻结 `process_input_data` 的角度/level 参数计算
  interface、normal 与 tangent bounds；没有把 bounds 改写成另一种优化定义。
- Regression smoothing 不是向对角线加 nugget。冻结源码先完成全部核/polynomial
  装配，再把 Single Surface 的 inequality+interface 对角项、Lajaunie 的 same-level
  increment 对角项直接替换为 `K((0,0,0),(0,0,smoothing_amount))`。即使 amount 为零
  也执行替换；Stratigraphic、Continuous Property 和 Vector Field 忽略该参数。
- `Continuous_Property` 的私有 polynomial helpers 在冻结可达参数路径中因
  `poly_term=false/n_poly_terms=0` 不被调用，且其 `.size()` 维数检查存在源级问题；T17
  只装配实际可达的 interface value block，没有借机臆增 polynomial 能力。其最终可达性
  仍由既有 T27 复核。
- `AssemblyError` 保留冻结 matrix/equality/inequality 阶段类别，并让 `R/AR` 缺失导数
  继续作为 `KernelError::LinearDerivativeUnavailable` 可见；modified-basis layout 与
  普通核（或反向）的不一致被安全拒绝，不让调用者得到静默错误矩阵。

T17 冻结 probe 对五模型普通路径、Single/Stratigraphic 普通 QP、三模型 restricted-range
和两种 smoothing 路径的完整矩阵与所有 RHS/bounds 做 row-major binary64 位哈希；Rust
全部精确匹配，并另行检查标签、分区、P/PT/0、符号、数值对称性和 smoothing 边界。
本任务不求解系统，不宣称 T18–T20 solver 或 T32 全局 parity 已完成。

## T18 partial-pivot LU 与后验判定

- `Linear_LU_decomposition::validate_matrix_systems` 只检查 interpolation matrix 是否全部
  有限；冻结模型路径不调用这个函数。`solve` 先检查 RHS 行数，再无条件调用
  `partialPivLu().solve`，且只因最终 weights 非有限而返回 `false`。debug-only 的 condition
  number 输出不参与接受/拒绝，`check_solution` 只打印 relative L2 residual 并恒定返回
  `true`。
- Rust 对动态小矩阵保留冻结 Eigen 的 unblocked partial-pivot 控制流：每列选择首个最大
  绝对值主元、交换整行、以 `a -= l*u` 次序更新尾块，并按 column-major dynamic vector
  路径依次更新 UnitLower/Upper RHS。输出保留 row transpositions、Eigen permutation indices、
  packed LU、pivot 值和首个 exact-zero pivot 作为可审查证据。
- 不以 condition number 或 pivot 大小预拒绝。3×3 Hilbert 系统实际进入求解并与冻结
  weights/pivots/packed LU 逐 bit 一致。冻结 Eigen 的 Upper solve 在当前 RHS 分量精确为零时
  跳过除法，因此 `[[1,2],[2,4]]/[3,6]` 虽有 exact-zero pivot，仍产生有限 `[3,0]`、零
  residual 并成功；Rust 忠实保留这一可观察的非唯一有限解语义。相同矩阵配不一致 RHS
  会在尝试后产生非有限 weights，并映射为 `SingularSystem` 与外层
  `Error::LinearSolverFailure`。
- 非有限 matrix/RHS 仍执行冻结式 factorization/solve 以保存 `attempted=true` 与 pivot
  evidence，随后分别稳定分类为 `NonFiniteMatrix`/`NonFiniteRightHandSide`；非法 storage、
  空系统、非方阵和 RHS 维数不符在 Rust 边界安全拒绝且不复制 Eigen release assertion/UB。
  `surfe_matrix_system_valid` 与完整 safe preflight 被分开暴露，避免把强化输入检查误写成
  冻结验证语义。
- 有限 weights 后计算 L2、relative L2、L-infinity residual 与 scale-aware backward error；
  接受门槛为 `64 * EPSILON * n`。这是求解后的 residual/constraint-feasibility 判定，不是
  条件数门槛；失败仍保留 candidate weights 和 residual evidence。良态、Hilbert、奇异一致
  与确定性对角占优系统均通过。T18 不实现 T19/T20 的 QP/LOQO。

## T19 普通 predictor-corrector QP 与后验可行性

- `Quadratic_Predictor_Corrector` 将 interpolation matrix 逐项乘二作为 Hessian，然后调用
  `Math_methods::quadratic_solver`；冻结模型路径注释掉 `validate_matrix_systems()` 调用。
  因此 Rust 单独保留 interpolation matrix 的有限性和 LLT 正定 evidence，但不把该结果
  当作求解门槛。冻结 indefinite 与 `1e-14` 病态样例都实际进入 KKT/LU 路径，Rust 保持
  相同分支和最终 binary64 weights。
- 初始化严格保留零 `x/y/z/s`、`sqrt(H.maxCoeff())` 的 `z/s`、首次完整 affine step，随后
  对 complementary variables 加 `1000 + 2*max_violation`；每轮 KKT、residual、predictor、
  affine step、`mu_aff`、立方 `sigma`、corrector 和最终 step 均保持冻结循环及乘加次序。
  `_find_step_length` 的符号分支、严格 `>1e-14` 和 unit cap 也逐控制流保留。LOQO 专用
  `_find_step/_find_positivity_step` 已核对边界，但实现仍严格留给既有 T20。
- 冻结停止顺序先检查 `iter > 5 && mu > prev_mu`，再检查 `mu < 1e-8`。Rust 保留两条分支，
  并提供每轮 `mu`、目标、stationarity、equality/inequality residual、实际最小 slack、两种
  step、`mu_aff` 和 `sigma`，以及每次 T18 KKT pivot/residual evidence。默认额外的 10000 轮
  上限只防止冻结无上限循环；达到上限是明确失败，不放宽停止规则制造成功。
- 冻结函数在 `nc=0` 时会除以零并索引空 `max_violation_list`；Rust 对缺少 inequality 的
  ordinary-QP 请求安全返回 `MissingInequalities`。其他空/维数/非有限错误也安全分类，
  外层稳定映射 `Error::PredictorCorrectorSolverFailure`，不复制 Eigen assertion 或无限循环。
- 冻结不可行用例 `x=0` 且 `x>=1` 在 `mu` 从 `1005418.6277759355` 变为
  `-1006.4220365718757` 后错误地返回成功、weights `[0]`、slack `-1`；该可观察行为先由
  oracle 固定。另一个 indefinite 用例在第六轮 `mu` 从 `17.105335153876929` 回升至
  `17.386342564804785` 后也返回不满足约束的成功 candidate。Rust 精确保留停止分支、trace
  和 candidate，再依据 T04 residual/feasibility 规则返回 `InfeasibleSolution`，不把源缺陷
  宣称为有效成功。
- 对有限且可行的非唯一系统不因退化预拒绝：冻结 zero-Hessian 样例产生 `[1,1000]`、零
  objective 和零终止 `mu`，Rust 逐 bit 保留并以后验 residual/feasibility 接受。普通 QP
  没有改用通用优化 crate，也没有实现 restricted-range/LOQO；后者仍仅归 T20。

## T20 restricted-range / LOQO 风格 QP

- `Quadratic_Predictor_Corrector_LOQO` 将 interpolation matrix 乘二为 `H`，直接求解
  `min 1/2*x'Hx`，约束保持冻结布局 `b <= A*x <= b+r`。因此对外报告的目标仍为
  `x' interpolation x`；`r` 是从下界到闭上界的非负宽度语义，不被改写成另一种约束、
  罚项或解后裁剪。T06 的 `SetRestrictedRange` 三参数原样写入 parameters，T17 的三模型
  bounded assembly 则提供这里直接消费的 `A/b/r`。
- 初始化保留冻结 `[-(H+I), A'; A, I]` KKT、右端 `[0;b]`，以及
  `g/z/t/s=max(abs(x),100)`、`v/w=max(abs(y),100)`、
  `p=max(abs(r-w),100)`、`q=v`。每轮继续使用 `D/E` 对角消元、predictor/corrector 两次
  T18 partial-pivot LU、四组 primal/dual positivity divisor 和固定 `0.95` fraction-to-boundary。
  中心参数使用 `((max(alpha_p,alpha_d)-1)/(max(alpha_p,alpha_d)+10))^2`；`P/Q` 两种
  对角比值的方向和所有更新次序均按冻结源码保留。
- 成功只沿 `significant_figures > 6` 分支；`dual_obj > primal_obj` 和 significant figures
  回退是冻结失败。源码循环无上限，Rust 默认增加 10000 轮安全上限且达到上限明确失败，
  不放宽成功规则。每轮公开目标、gap-derived significant figures、滞后一轮的 primal/dual
  infeasibility、predictor/corrector divisor、fraction 与 `mu`；每次 KKT 公开 stage、pivot
  transpositions 和 residual，终点另公开约束残差、最小上下界 slack、互补性和接受阈值。
- 冻结 `validate_matrix_systems()` 只记录 `H` 有限且 LLT 正定，但模型 `solve()` 不调用它；
  Rust 同样把此结果作为 evidence，而不是 condition-number/positive-definite gate。`1e-14`
  病态和 zero-Hessian 用例都会实际求解；zero-Hessian 的有限负极小 candidate 与负零目标
  逐 bit 保留。非有限初始 KKT 的 candidate 也保留到第一次 predictor 检查，再安全映射为
  `NonFiniteInput`，而不复制后续 Eigen 非有限运算。
- 冻结 tight zero-width、negative-range 和测试中的 indefinite 输入分别沿原始
  `dual_obj > primal_obj` 失败；Rust 保存相同 trace/候选/残差并映射稳定
  `Error::LoqoSolverFailure`。非法空/非方阵/维数在 Rust 边界安全拒绝。成功候选还必须通过
  有限性、stationarity/约束残差、闭区间上下界和互补性后验检查；这些检查不会用条件数
  提前改变冻结迭代分支。

T20 冻结 probe 覆盖 inactive、lower-active、两变量双边、tight、病态、zero-Hessian、
indefinite、negative-range 和 non-finite 九类；Rust 的成功 weights/objectives、停止迭代、
关键 `mu`/目标/失败分支与该 oracle 一致，并增加解析 upper-active 镜像检查。T20 只交付
QP 求解与 `A/b/r` 语义，不执行 T21 的活跃约束选择、普通 RBF 重装配或二次 LU。

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
