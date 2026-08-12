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
