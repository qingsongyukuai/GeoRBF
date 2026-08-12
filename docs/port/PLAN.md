# GeoRBF：Surfe 纯 Rust 忠实移植固定计划

## 目标与约束

源仓库固定为 `https://github.com/MichaelHillier/surfe.git`，参考提交固定为 `290dbe0ab344f4258a4935f05cad0f153f0f69a4`。本计划仅迁移 Surfe 非可视化核心；排除 Qt、VTK、GUI、`geo_builder`、等值面显示和纯可视化代码。

最终 GeoRBF 生产核心必须为纯 Rust。正常 Cargo 构建、测试和发布不得依赖 C++、Eigen、Qt、VTK、CMake、bindgen、CXX bridge、BLAS/LAPACK/MKL/OpenBLAS FFI。冻结 C++ Surfe 只允许作为仓库外部或被忽略目录中的 reference oracle，不得进入生产源码和发布包。

这是固定任务序列，不是动态 Issue 生成器。每个会话只处理 `STATE.json` 指定的一个当前任务；当前任务未满足门槛时不得前进。不得创建新的顶级任务、Issue、子 Issue、PR、分支或额外 worktree，不得改变顺序或扩大范围，除非用户明确要求。每个完成任务必须连同状态和日志形成一个原子提交。

下列 GeoRBF 模块路径是计划目标；若前序仓库审查发现已有兼容模块，应优先复用并在对应任务日志记录最终路径，不得仅因计划名称不同而重写。

## 固定任务序列

### 基线与判定体系

## T00 — 冻结 Surfe reference、许可证、迁移范围、兼容策略和纯 Rust 边界

- 任务编号：T00。
- 目标：固定源码身份、许可证义务、纳入/排除边界、兼容原则、reference 搜索优先级与纯 Rust 发布边界，为后续任务提供不可漂移的基线。
- 依赖：INIT。
- 对应 Surfe 文件/函数：`License.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`README.md@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`.gitmodules@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`CMakeLists.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；入口 `Surfe_API::ComputeInterpolant`、`EvaluateInterpolantAtPoint(s)`、`EvaluateVectorInterpolantAtPoint(s)`。
- GeoRBF 目标模块：`docs/port/source-baseline.md`、`docs/port/compatibility.md`、许可证与依赖边界文档。
- 具体交付物：固定提交校验记录；许可证和来源追踪规则；纳入/排除清单；有效输入、缺陷兼容和错误分类策略；reference 路径解析协议；生产/测试/发布的纯 Rust 禁止依赖清单。
- 必须运行的验证：`git -C "$SURFE_REFERENCE_DIR" rev-parse HEAD` 或按规定优先级解析出的等价命令；许可证文本校验；`git diff --check`；检查文档同时包含固定提交、排除项、reference 优先级和全部 native 禁止项。
- 完成门槛：reference 可获得且提交精确匹配；许可证与边界文档无歧义；所有纳入和排除项均可判定；未把源码或构建产物加入仓库。
- 禁止事项：不得复制或提交 Surfe 源码；不得建立 oracle、golden fixture 或实现算法；不得把 GUI/VTK/Qt/`geo_builder` 纳入迁移；reference 不可用时不得伪造通过。

## T01 — 建立完整 C++ 清单、调用链、数据流和依赖图

- 任务编号：T01。
- 目标：枚举迁移范围内全部 C++ 文件、类型、函数、调用关系、数据流和第三方依赖，形成可审计覆盖基线。
- 依赖：T00。
- 对应 Surfe 文件/函数：`surfe_lib/*.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`math_lib/*.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_pybindings/*@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 中到核心的适配调用、`test/main.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；重点入口 `Surfe_API`、`GRBF_Modelling_Methods::get_method`、`setup_basis_functions`、各模型 `process_input_data/get_interpolation_matrix/setup_system_solver/eval_*`、`System_Solver::solve`。
- GeoRBF 目标模块：`docs/port/inventory.md`、`docs/port/call-graph.md`、`docs/port/data-flow.md`。
- 具体交付物：逐文件和逐符号清单；五种模型从 API 到求解与评估的调用链；约束到矩阵、RHS、权重和预测的字段级数据流；Eigen/CMake/可视化依赖边界图；未决项只能指向已有任务。
- 必须运行的验证：用 `rg`/符号提取命令将冻结源码声明与清单互相核对；检查每个迁移范围内 `.h/.cpp` 至少出现一次；检查所有工厂分支和虚函数覆写均有归属；`git diff --check`。
- 完成门槛：范围内文件、类型和函数覆盖无缺口；入口、调用者、被调用者和数据所有权清楚；所有依赖均标为迁移、Rust 替代或排除。
- 禁止事项：不得实现任何 Rust 算法；不得把声明、TODO 或不可达函数当成可用能力；不得创建动态任务或改变固定序列。

## T02 — 逐符号能力分类

- 任务编号：T02。
- 目标：依据定义、调用点和可达性证据，把 T01 的每个符号分类为已实现、部分实现、TODO、不可达、缺陷、可视化或排除。
- 依赖：T01。
- 对应 Surfe 文件/函数：T01 清单中的全部符号，特别是 `basis.cpp` 的 `WendlandC2/MaternC4/MQ3`，`continuous_property.cpp` 的方法覆写，`modeling_methods.cpp` 的 `run_greedy_algorithm`，各模型 `convert_modified_kernel_to_rbf_kernel`，以及 `geo_builder/*@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的排除边界。
- GeoRBF 目标模块：`docs/port/symbol-classification.md`。
- 具体交付物：逐符号唯一分类表；定义与调用证据；部分实现的缺口；TODO/不可达/缺陷的可观察含义；可视化/排除理由；每项映射到一个且仅一个既有任务或明确排除。
- 必须运行的验证：以冻结源码声明集和定义集做双向集合核对；搜索 `TODO`、空函数、恒定返回、抛异常和未调用符号；用调用链验证可达性；`git diff --check`。
- 完成门槛：T01 的每个范围内符号恰有一个有证据的分类；不存在“未分类”或仅凭名称推断；不可达与缺陷没有被误报为能力。
- 禁止事项：不得实现分类结果；不得新增任务；不得把测试未覆盖解释为源码不可达；不得因源码缺陷静默定义新语义。

## T03 — 建立仓库外冻结 Surfe reference oracle

- 任务编号：T03。
- 目标：在仓库外部或被忽略目录构建可重复调用的冻结 C++ oracle，仅向仓库提交协议、字段清单和生成说明。
- 依赖：T02。
- 对应 Surfe 文件/函数：`surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的约束设置、`ComputeInterpolant`、标量/向量评估；`surfe_lib/{basis,matrix_solver,modeling_methods,single_surface,lajaunie,stratigraphic_surfaces,continuous_property,vector_field}.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的中间值提取点；`test/main.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`。
- GeoRBF 目标模块：`docs/port/oracle-protocol.md`、`docs/port/oracle-manifest.json`、只含协议验证的 Rust 测试辅助模块。
- 具体交付物：reference 发现与提交校验说明；确定性请求/响应协议；核、导数、矩阵、RHS、求解证据、预测和错误的字段清单；构建与运行说明；源码和生成物的忽略规则。
- 必须运行的验证：在解析出的冻结 reference 上构建并运行最小 oracle smoke case；校验响应 schema、确定性和 commit 身份；`git status --ignored` 确认 C++ 与二进制未被跟踪；`git diff --check`。
- 完成门槛：同一输入重复输出稳定；协议能覆盖后续 parity 所需中间值和错误；仓库未跟踪任何 Surfe 源码、对象或二进制；无法运行的 oracle 不得标为完成。
- 禁止事项：不得提交 C++、Eigen、CMake 产物或正式 golden fixture；不得让 Cargo 的正常构建/测试依赖 oracle；不得修改冻结 Surfe 算法以制造一致。

## T04 — 定义 golden fixture 与差分测试协议

- 任务编号：T04。
- 目标：定义确定性、可追踪、可版本化的 golden fixture 格式、覆盖数据集、分层容差和差分判定协议。
- 依赖：T03。
- 对应 Surfe 文件/函数：oracle 暴露的 `RBFKernel::{basis,dx_p1,dx_p2,dy_p1,dy_p2,dz_p1,dz_p2,dxx..dzz}`、`Modified_Kernel::basis_*`、各模型 `get_interpolation_matrix/get_equality_values/get_inequality_*`、求解器 `solve`、`eval_scalar_interpolant_at_point`、`eval_vector_interpolant_at_point`。
- GeoRBF 目标模块：`docs/port/fixtures.md`、`tests/fixtures/schema/`、`tests/common/parity` 的格式读取与比较基础设施。
- 具体交付物：带 source commit、平台、精度和生成命令的 fixture schema；确定性数据集清单；离散字段精确比较规则；数值量分层容差；NaN/Inf/错误编码规则；权重非唯一时的残差/可行性/预测判定；再生成和审阅流程。
- 必须运行的验证：schema 正反例测试；序列化往返测试；同一 oracle 输入重复生成字节稳定性测试；容差边界单元测试；`cargo test`（若 Cargo 工程已存在）和 `git diff --check`。
- 完成门槛：协议覆盖核、两点一阶导数、混合 Hessian、anisotropy、Modified Kernel、矩阵/RHS、LU/QP、标量场、梯度场、目标值和迭代证据；离散字段无模糊容差；本任务尚不生成正式数值 fixture。
- 禁止事项：不得实现迁移算法；不得提交未审阅的正式 golden 数据；不得用单一宽松容差掩盖符号、布局或求解分支差异。

### Rust 数学基础

## T05 — 审查 GeoRBF 并建立纯 Rust 工程护栏

- 任务编号：T05。
- 目标：先审查和复用现有 GeoRBF，再建立纯 Rust 依赖、CI、格式、静态检查和 native dependency 审计护栏。
- 依赖：T04。
- 对应 Surfe 文件/函数：`CMakeLists.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`.gitmodules@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/CMakeLists.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`（若存在于清单）、Eigen 类型使用点及 `surfe_lib_module.h` 导出边界；不迁移构建系统本身。
- GeoRBF 目标模块：现有 Cargo workspace、`.github/workflows/`、依赖审计脚本/测试、`docs/port/georbf-audit.md`。
- 具体交付物：现有模块和可复用性审计；目标模块图；最小纯 Rust 依赖选择及理由；fmt/clippy/test/doc/多平台 CI；禁止 native 依赖和 build script 的机器可检护栏；发布包内容审计基线。
- 必须运行的验证：`cargo fmt --check`、`cargo clippy --all-targets --all-features -- -D warnings`、`cargo test --all-targets --all-features`、`cargo tree --all-features`、`cargo package --list`，以及针对 `build.rs`、`links`、`cc`、`cmake`、`bindgen`、`cxx`、BLAS/LAPACK/MKL/OpenBLAS 的审计命令。
- 完成门槛：现有代码逐项判定复用/扩展/替换且替换有证据；CI 与本地护栏能拒绝所有禁止依赖；所有验证通过；未引入算法实现。
- 禁止事项：不得无理由重建工程或公共 API；不得添加 C/C++ 构建桥；不得实现 T06 及以后数学功能。

## T06 — 参数、枚举、默认值、名称和类型化错误

- 任务编号：T06。
- 目标：忠实迁移 Surfe 参数模型、枚举判别、默认值、核/模型名称解析和可观察错误分类。
- 依赖：T05。
- 对应 Surfe 文件/函数：`surfe_lib/modelling_parameters.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Parameter_Types::{DWRT,SecondDerivatives,FirstDerivatives,RBF,SolverType,ModelType,AXIS}`、`Parameters/InternalParameters/InputParameters`；`surfe_lib/grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 全部异常；`surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `SetRBFKernel`、参数 setter、`get_method_from_parameters`。
- GeoRBF 目标模块：`src/parameters.rs`、`src/error.rs`、名称兼容测试。
- 具体交付物：Rust 枚举和参数结构；逐字段默认值；Surfe 兼容字符串与拒绝规则；稳定的类型化错误类别及源异常映射；不泄漏 C++ 表示的安全构造接口。
- 必须运行的验证：参数默认值 golden comparison；所有合法/非法名称表驱动测试；错误分类测试；Serde/显示往返测试（若暴露）；`cargo fmt --check`、`cargo clippy --all-targets --all-features -- -D warnings`、`cargo test --all-targets --all-features`。
- 完成门槛：离散默认、枚举、名称和接受/拒绝行为精确匹配；每个范围内 Surfe 异常有明确 Rust 类别；无算法或矩阵行为提前实现。
- 禁止事项：不得为“更友好”接受 Surfe 不接受的名称；不得丢失错误类别；不得实现约束、核或求解器。

## T07 — 约束数据类型与方向转换

- 任务编号：T07。
- 目标：移植 `Point`、`Interface`、`Planar`、`Tangent`、`Inequality` 的数据语义及 strike/dip/azimuth/polarity/normal 方向转换。
- 依赖：T06。
- 对应 Surfe 文件/函数：`surfe_lib/modelling_input.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Point/Interface/Inequality/Planar/Tangent`，`Planar::_compute_strike_dip_polarity_from_normal`、`_compute_normal_from_strike_dip_polarity`、`getDipVector`、`getStrikeVector`、`setNormalBounds`。
- GeoRBF 目标模块：`src/geometry.rs`、`src/constraints.rs`。
- 具体交付物：安全 Rust 值类型；构造与访问语义；方向转换和角度单位；法向量/切向量校验；边界计算字段；针对零向量和非有限输入的明确错误。
- 必须运行的验证：冻结 C++ 的方向转换 golden；往返与单位向量解析恒等式；极角、极性、象限和边界用例；有限差分不适用但须做几何恒等式；标准 fmt/clippy/test 全套。
- 完成门槛：有效输入转换数值在既定容差内等价；离散极性和角度约定精确；安全拒绝不复制 UB；公共 API 尚未扩张到 fitted model。
- 禁止事项：不得擅自改用不同地质方位约定；不得以归一化“修复”改变 Surfe 有效输入结果；不得实现排序或空间算法。

## T08 — 排序、同点判定、去重与 level/reference 分组

- 任务编号：T08。
- 目标：忠实移植排序、`1e-3` 同点判定、分类去重、精确 level 分组和 interface reference point 分组。
- 依赖：T07。
- 对应 Surfe 文件/函数：`surfe_lib/modelling_input.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Point::operator<`、`operator==`；`surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `_get_distinct_interface_iso_values`、`_get_interface_points`、`_get_distinct_inequality_iso_values`、`get_interface_data`、`remove_collocated_constraints`；`math_lib/math_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Math_methods::sort_vector_w_index`。
- GeoRBF 目标模块：`src/ordering.rs`、`src/constraints/grouping.rs`。
- 具体交付物：确定性总排序适配；`1e-3` 同点谓词；按约束类别的保留/移除规则；精确 level 分组；reference point 选择和分组索引；重复和边界数据集。
- 必须运行的验证：排序/索引 golden；恰在 `1e-3` 两侧的同点测试；跨类别和同类别重复矩阵；`-0.0`、相等 level、输入排列测试；分组离散结果精确比较；标准 fmt/clippy/test 全套。
- 完成门槛：排序、去重、level 组、reference point 及自由度前置计数精确一致并确定性；没有近似 level 合并。
- 禁止事项：不得用哈希迭代顺序影响输出；不得将 level 改为容差分组；不得实现增量 pair 或矩阵装配。

## T09 — 空间辅助算法

- 任务编号：T09。
- 目标：移植距离、最近邻、平均最近邻、bounds、极值点和相关空间辅助算法。
- 依赖：T08。
- 对应 Surfe 文件/函数：`surfe_lib/modelling_input.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `convert_constraints_to_points`、`distance_btw_pts`、`nearest_neighbour_index`、`get_n_nearest_neighbours_to_point`、两个 `furtherest_neighbour_index`、`avg_nn_distance`、`spatial_metrics`、`Find_STL_Vector_Indices_FurtherestTwoPoints`、`Find_STL_Vector_Index_ofPointClosestToOtherPointWithinDistance`、`calculate_bounds`、`get_extremal_point_data_indices_from_points`、`get_largest_distance_between_points`、`get_maximal_axial_variability_order` 及 `Constraints::compute_*_avg_nn_distance`。
- GeoRBF 目标模块：`src/spatial.rs`。
- 具体交付物：纯 Rust 空间函数；空集/单点/并列距离语义；确定性索引 tie-break；`SpatialParameters` 等价结果；残差候选筛选辅助函数的归属记录。
- 必须运行的验证：C++ golden；手算距离/bounds；平移和旋转相关性质；并列、重复、退化轴、空/单点错误测试；标准 fmt/clippy/test 全套。
- 完成门槛：所有离散索引精确一致，数值空间量在容差内一致；边界错误有类型化结果且不复制越界 UB。
- 禁止事项：不得引入 k-d tree 或近似最近邻改变 tie-break；不得提前实现 Greedy 策略。

## T10 — Polynomial 与 truncated polynomial

- 任务编号：T10。
- 目标：移植 0/1/2 阶 polynomial 以及模型所用 truncated polynomial 的值和一阶导数，固定项顺序。
- 依赖：T09。
- 对应 Surfe 文件/函数：`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Polynomial_Basis`、`Poly_Zero/Poly_First/Poly_Second::{basis,dx,dy,dz}`；各模型 `_get_polynomial_matrix_block` 中的截断使用和项截取。
- GeoRBF 目标模块：`src/polynomial.rs`。
- 具体交付物：固定长度或安全向量表示；0/1/2 阶项序；值与 x/y/z 导数；truncation 规则；模型矩阵尚未装配的独立测试接口。
- 必须运行的验证：C++ golden；解析导数；有限差分；项顺序和维数精确测试；多尺度/负坐标/原点测试；标准 fmt/clippy/test 全套。
- 完成门槛：值、导数、项序、维数和截断精确满足三角验证；零值和常数项行为清楚。
- 禁止事项：不得重排为更常见的 monomial 顺序；不得引入高阶项；不得装配模型矩阵。

## T11 — Unisolvent 选择与 Lagrangian polynomial basis

- 任务编号：T11。
- 目标：移植四点 unisolvent 选择、退化处理和 Lagrangian polynomial basis。
- 依赖：T10。
- 对应 Surfe 文件/函数：`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Lagrangian_Polynomial_Basis::_get_unisolvent_subset`、`_initialize_basis`、`poly`、`poly_dx`、`poly_dy`、`poly_dz`；`failurecreatinglagrangianpolynomialbasis`。
- GeoRBF 目标模块：`src/polynomial/lagrangian.rs`。
- 具体交付物：确定性四点选择；退化和不足点错误；basis 系数与评估；导数；所选原始索引证据；与 Modified Kernel 解耦的测试。
- 必须运行的验证：C++ 选择索引与数值 golden；Kronecker 插值恒等式；解析/有限差分导数；共面、共线、重复、近退化和输入排列测试；标准 fmt/clippy/test 全套。
- 完成门槛：对有效输入选择与 basis 等价；退化失败类别等价且无未定义内存访问；所有三角验证通过。
- 禁止事项：不得用随机选点；不得以不同 pivot 规则改变有效输入选择；不得实现 Modified Kernel。

## T12 — 九种 isotropic 核及导数

- 任务编号：T12。
- 目标：移植冻结源码定义的九种 isotropic 核、对第一/第二点的一阶导数、混合二阶导数和零距离行为。
- 依赖：T11。
- 对应 Surfe 文件/函数：`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `RBFKernel::{radius,basis_pt_pt,basis_pt_planar_*,basis_planar_*_pt,basis_planar_planar}` 及 `Cubic`、`Gaussian`、`MQ`、`MQ3`、`TPS`、`IMQ`、`R`、`WendlandC2`、`MaternC4` 的 `basis`、`dx_p1..dz_p2`、`dxx..dzz`。
- GeoRBF 目标模块：`src/kernel/isotropic.rs`、`src/kernel/derivatives.rs`。
- 具体交付物：九种核的纯 Rust 实现；shape/support 参数；两点导数方向；完整混合 Hessian；零距离/支撑边界分支；统一内部核接口。
- 必须运行的验证：每个核的冻结 C++ golden；径向、交换、符号和 Hessian 对称解析恒等式；对两个点分别做有限差分；零距离、近零、支撑内外和多尺度测试；标准 fmt/clippy/test 全套。
- 完成门槛：九种核全部通过三角验证；离散分支精确一致；任何源码/数学冲突已记录且标准有效输入保持 parity。
- 禁止事项：不得只用自动微分替代源码语义而不做 parity；不得更改零距离分支；不得实现 anisotropy 或 Modified Kernel。

## T13 — 全局各向异性

- 任务编号：T13。
- 目标：移植 Surfe 的全局各向异性计算、`f32` 中间语义、特征值截断和支持矩阵。
- 依赖：T12。
- 对应 Surfe 文件/函数：`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `RBFKernel::get_global_anisotropy`、`scaled_radius`、`_Global_Plunge`、Eigen `Matrix3f _Transform`，以及 `ACubic/AGaussian/AMQ/ATPS/AIMQ/AR` 的值与导数路径。
- GeoRBF 目标模块：`src/kernel/anisotropy.rs`。
- 具体交付物：纯 Rust 3×3 支持矩阵/特征分解路径；明确 `f32` 中间舍入；特征值排序、截断和失败规则；变换后的半径、核值和导数。
- 必须运行的验证：变换、特征值和 plunge 的 C++ golden；正交/对称/尺度解析检查；各向同性极限；有限差分核导数；退化法向和特征值边界；标准 fmt/clippy/test 全套。
- 完成门槛：`f32` 中间导致的可观察结果在分层容差内等价；矩阵布局和特征值截断分支一致；失败类别明确。
- 禁止事项：不得直接改成全 `f64` 并声称等价；不得借助 native BLAS/LAPACK；不得提前实现 Modified Kernel。

## T14 — Modified Kernel 全组合

- 任务编号：T14。
- 目标：移植 Modified Kernel 的值、梯度、Hessian、Planar 和 Tangent 全部组合。
- 依赖：T13。
- 对应 Surfe 文件/函数：`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Modified_Kernel` 构造及 `basis_pt_pt`、`basis_pt_planar_*`、`basis_planar_*_pt`、`basis_pt_tangent`、`basis_tangent_pt`、`basis_planar_planar`、`basis_tangent_tangent`、`basis_planar_tangent`、`basis_tangent_planar`。
- GeoRBF 目标模块：`src/kernel/modified.rs`。
- 具体交付物：以 T11 Lagrangian basis 和 T12/T13 核为基础的纯 Rust Modified Kernel；Value/Planar/Tangent 全笛卡尔组合；梯度/Hessian 中间证据；构造失败映射。
- 必须运行的验证：所有组合 C++ golden；消去/正交解析性质；参数交换和混合导数对称性；有限差分；isotropic/anisotropic 和零距离矩阵；标准 fmt/clippy/test 全套。
- 完成门槛：值、梯度、Hessian、Planar、Tangent 全组合均通过三角验证；无漏分支；错误分类等价。
- 禁止事项：不得只实现模型当前首个调用到的组合；不得改变 Lagrangian 选择；不得装配完整系统。

### 泛函、矩阵和求解器

## T15 — 统一内部线性泛函

- 任务编号：T15。
- 目标：建立 Value、Derivative、Tangent、Difference 的统一内部线性泛函，精确表达冻结 Surfe 的符号和对两个核参数的作用方向。
- 依赖：T14。
- 对应 Surfe 文件/函数：`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的全部 `Kernel::basis_*`；`lajaunie.cpp`、`stratigraphic_surfaces.cpp` 的 increment pair 差值块；各模型 `get_interpolation_matrix` 中的 point/planar/tangent 组合。
- GeoRBF 目标模块：`src/functional.rs`。
- 具体交付物：带作用点和方向的 `Value/Derivative/Tangent/Difference` 表示；核双线性作用 API；差值展开；自由度标签；普通核与 Modified Kernel 的统一调用层。
- 必须运行的验证：所有泛函对 C++ `basis_*` golden；线性、差值展开、参数交换和符号恒等式；有限差分校验 derivative/tangent；标准 fmt/clippy/test 全套。
- 完成门槛：每个 Surfe 核组合有唯一无歧义表示；符号和作用点精确；不依赖模型特定索引即可测试。
- 禁止事项：不得在泛函层重排约束；不得合并数学上相似但 Surfe 符号不同的分支；不得装配模型布局。

## T16 — 五种模型的确定性 constraint layout

- 任务编号：T16。
- 目标：建立五种模型各自确定性的 constraint layout、索引顺序、自由度计数和矩阵分区。
- 依赖：T15。
- 对应 Surfe 文件/函数：`single_surface.cpp`、`lajaunie.cpp`、`stratigraphic_surfaces.cpp`、`continuous_property.cpp`、`vector_field.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `process_input_data`、`get_method_parameters`、`get_interpolation_matrix`、`get_equality_values`、`get_inequality_*`；`modeling_methods.h` 内部参数字段。
- GeoRBF 目标模块：`src/layout.rs`、`src/model/*/layout.rs`。
- 具体交付物：五种 layout 描述；每个区段的约束来源、行列范围和自由度；polynomial 分区；equality/inequality 分区；稳定的索引查询与调试快照。
- 必须运行的验证：对各类约束组合做 C++ 维数/行列标签 golden；空类别和多 level 情况；输入排列和去重后的 layout；精确快照测试；标准 fmt/clippy/test 全套。
- 完成门槛：五种模型的维数、顺序、分区和自由度计数精确一致且确定性；尚不填充数值矩阵。
- 禁止事项：不得为统一性改变模型原顺序；不得使用 unordered/hash 顺序；不得装配数值或调用求解器。

## T17 — 统一矩阵、RHS 与 smoothing 装配

- 任务编号：T17。
- 目标：实现普通核块、polynomial block、modified block、RHS 和 regression smoothing 的统一装配，并保持各模型布局与符号。
- 依赖：T16。
- 对应 Surfe 文件/函数：五个模型的 `get_interpolation_matrix`、`_get_polynomial_matrix_block`、`_insert_polynomial_matrix_blocks_in_interpolation_matrix`、`get_equality_values`、`get_inequality_matrix`、`get_inequality_values`；`modeling_methods.cpp::get_equality_matrix`；参数 `use_regression_smoothing/smoothing_amount`。
- GeoRBF 目标模块：`src/assembly.rs`、`src/model/*/assembly.rs`。
- 具体交付物：稠密纯 Rust 矩阵/向量表示；按 T16 layout 装配所有核/多项式/modified 块；RHS；smoothing 对角或源码指定项；中间矩阵快照接口。
- 必须运行的验证：逐块和完整矩阵/RHS C++ golden；维数、行列标签、对称性和符号检查；无 smoothing/边界值测试；五模型最小数据集；标准 fmt/clippy/test 全套。
- 完成门槛：完整矩阵与 RHS 在分层容差内 parity，离散布局精确；所有模型所需块均覆盖；不求解系统。
- 禁止事项：不得改用稀疏/低秩/局部近似；不得因条件数拒绝装配；不得实现 LU/QP。

## T18 — Pure-Rust partial-pivot LU

- 任务编号：T18。
- 目标：实现纯 Rust partial-pivot LU 路径及与 Surfe 一致的验证、尝试求解和失败语义。
- 依赖：T17。
- 对应 Surfe 文件/函数：`surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `System_Solver`、`Linear_LU_decomposition::{solve,validate_matrix_systems,check_solution}`，及 Eigen `partialPivLu().solve` 的可观察路径。
- GeoRBF 目标模块：`src/solver/lu.rs`、`src/solver/error.rs`。
- 具体交付物：纯 Rust partial-pivot LU 分解/求解；矩阵/RHS 验证；有限性、残差和可行性后验检查；失败分类；pivot 与残差证据。
- 必须运行的验证：冻结 Surfe 的良态、病态、奇异和非有限系统 golden；已知解与随机确定性系统；残差和 backward-error 检查；确认病态系统先尝试求解；标准 fmt/clippy/test 和 native 依赖审计。
- 完成门槛：正常路径预测/残差 parity；病态不因条件数预拒绝；奇异/非有限失败分类符合兼容策略；生产依赖纯 Rust。
- 禁止事项：不得调用 LAPACK/BLAS/MKL/OpenBLAS FFI；不得以条件数阈值代替 Surfe 求解语义；不得实现 QP。

## T19 — 普通 predictor-corrector QP

- 任务编号：T19。
- 目标：忠实移植普通 predictor-corrector QP 的初始化、KKT 系统、方向、步长、停止和失败规则。
- 依赖：T18。
- 对应 Surfe 文件/函数：`surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Quadratic_Predictor_Corrector::{solve,validate_matrix_systems}`；`math_lib/math_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Math_methods::quadratic_solver`、`_find_step`、`_find_positivity_step`、`_find_step_length`、`max_element_wrt_zero`。
- GeoRBF 目标模块：`src/solver/qp_predictor_corrector.rs`。
- 具体交付物：纯 Rust primal-dual 状态；KKT 装配与 LU 调用；预测/校正方向；步长和中心参数；停止/迭代上限；目标值、残差、互补性和迭代 trace。
- 必须运行的验证：oracle 的逐迭代或关键迭代证据 golden；解析小 QP；等式+不等式、边界活跃、退化和不可行用例；KKT 残差/可行性/目标；标准 fmt/clippy/test 与 native 审计。
- 完成门槛：初始化和分支精确，收敛用例的目标/可行性/预测 parity；失败类别一致；权重可不逐位相同但残差、约束和预测必须满足全局规则。
- 禁止事项：不得换用通用优化 crate 的不同算法；不得放宽停止规则制造成功；不得实现 restricted-range/LOQO。

## T20 — Restricted-range / LOQO 风格 QP

- 任务编号：T20。
- 目标：忠实移植 restricted-range/LOQO 风格 QP 和上下界语义。
- 依赖：T19。
- 对应 Surfe 文件/函数：`surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Quadratic_Predictor_Corrector_LOQO::{solve,validate_matrix_systems}`；`math_lib/Math_methods::quadratic_solver_loqo` 及步长 helpers；`Surfe_API::SetRestrictedRange`。
- GeoRBF 目标模块：`src/solver/qp_loqo.rs`。
- 具体交付物：上下界/松弛变量语义；LOQO 初始化、KKT、步长、停止、失败；restricted range 的界映射；完整迭代证据。
- 必须运行的验证：C++ 逐阶段 golden；解析 box-constrained QP；单/双边界、紧边界、不可行、退化和病态用例；目标、残差、可行性、互补性；标准 fmt/clippy/test 与 native 审计。
- 完成门槛：界的开闭、符号和布局精确；有效输入的目标/预测/可行性 parity；失败分类一致且不因条件数预拒绝。
- 禁止事项：不得把 bounds 转成不同优化定义；不得换算法；不得提前做模型重建。

## T21 — QP 到普通 RBF 线性系统的重建

- 任务编号：T21。
- 目标：移植 QP 结果转普通 RBF 线性系统的重建路径，保留被选约束、目标值/level 更新和再次求解语义。
- 依赖：T20。
- 对应 Surfe 文件/函数：`single_surface.cpp`、`lajaunie.cpp`、`stratigraphic_surfaces.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `convert_modified_kernel_to_rbf_kernel`、`setup_system_solver`、`get_inequality_values`，以及相关 `process_input_data/get_interpolation_matrix` 重入路径。
- GeoRBF 目标模块：`src/model/reconstruct.rs`。
- 具体交付物：QP 解到等式约束/iso-value 的确定性转换；普通核 layout 重建；二次装配与 LU 求解；转换前后索引映射和预测证据。
- 必须运行的验证：C++ 的选中约束、重建矩阵/RHS、LU 残差和预测 golden；无活跃/多活跃边界；重建前后可行性；标准 fmt/clippy/test 全套。
- 完成门槛：重建的离散选择、布局和分支精确；最终场与 oracle parity；所有失败能区分 QP、重建装配和 LU。
- 禁止事项：不得直接复用 QP 权重跳过源码重建；不得改变活跃集定义；不得提前完成具体模型 API。

### 五种建模方法

## T22 — Single Surface 线性纵向切片

- 任务编号：T22。
- 目标：完成 Single Surface 的普通线性 equality 路径和纵向切片语义，包括约束处理、装配、LU 与标量/梯度评估。
- 依赖：T21。
- 对应 Surfe 文件/函数：`surfe_lib/single_surface.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的构造、`get_method_parameters`、`process_input_data`、`get_interpolation_matrix`、`get_equality_values`、`setup_system_solver`、`eval_scalar_interpolant_at_point`、`eval_vector_interpolant_at_point`；基类 `check_input_data/check_interpolant`。
- GeoRBF 目标模块：`src/model/single_surface.rs`。
- 具体交付物：Single Surface 普通构建/拟合内部路径；接口、Planar、Tangent 的线性系统；polynomial 与 smoothing；单点/批量内部评估；中间布局和残差证据。
- 必须运行的验证：oracle 的清洗后约束、layout、矩阵、RHS、权重残差、标量和梯度 golden；每种核/多项式代表用例；纵向切片数据；错误路径；标准 fmt/clippy/test 全套。
- 完成门槛：普通 equality 分支的离散行为精确、数值场 parity；不要求病态权重逐位相同，但残差和预测符合全局规则。
- 禁止事项：不得实现 inequality/QP、restricted range 或 Modified Kernel 分支；不得暴露 T29 最终公共 Builder。

## T23 — Single Surface inequality/QP

- 任务编号：T23。
- 目标：完成 Single Surface 的 inequality 约束和普通 predictor-corrector QP 路径。
- 依赖：T22。
- 对应 Surfe 文件/函数：`single_surface.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `get_inequality_matrix`、两个 `get_inequality_values`、`setup_system_solver`、`process_input_data`、`get_interpolation_matrix`，以及 `Quadratic_Predictor_Corrector` 调用链。
- GeoRBF 目标模块：`src/model/single_surface.rs` 的 inequality/QP 路径。
- 具体交付物：inequality 符号和 RHS；QP 分区；求解与结果验证；活跃约束证据；标量/梯度输出；错误映射。
- 必须运行的验证：上下 level inequality、混合 equality/inequality、活跃/非活跃和不可行数据的 C++ golden；矩阵/RHS、QP 目标、可行性和场；标准 fmt/clippy/test 全套。
- 完成门槛：不等式接受/拒绝、符号、布局、求解分支和预测 parity；失败分类准确。
- 禁止事项：不得实现 restricted-range/LOQO 或 Modified Kernel；不得用后处理裁剪替代 QP。

## T24 — Single Surface restricted-range、Modified Kernel 与线性重建

- 任务编号：T24。
- 目标：完成 Single Surface restricted-range、Modified Kernel、LOQO 和转普通 RBF 线性系统的完整路径。
- 依赖：T23。
- 对应 Surfe 文件/函数：`single_surface.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `setup_basis_functions` 相关分支、`get_inequality_values(VectorXd&,VectorXd&)`、`setup_system_solver`、`convert_modified_kernel_to_rbf_kernel`、评估函数；`Modified_Kernel` 与 `Quadratic_Predictor_Corrector_LOQO`。
- GeoRBF 目标模块：`src/model/single_surface.rs` 的 restricted/modified/reconstruct 路径。
- 具体交付物：range bounds；Modified Kernel 构建；LOQO 调用；活跃结果转换；普通核再装配/LU；最终 scalar/gradient；全路径 trace。
- 必须运行的验证：restricted range 边界、多个 interface level、anisotropy/普通核代表数据的端到端 C++ golden；modified 与重建中间矩阵；QP/LU 残差和最终场；标准 fmt/clippy/test 全套。
- 完成门槛：Single Surface 所有纳入分支可达且 parity；range 与重建离散选择一致；失败阶段可区分。
- 禁止事项：不得跳过重建；不得把 range 解释成输出裁剪；不得开始 Lajaunie。

## T25 — Lajaunie reference point、increment 与 iso-value

- 任务编号：T25。
- 目标：完成 Lajaunie 方法的 reference point、同层 increment pair、系统装配、求解、评估和 iso-value 更新。
- 依赖：T24。
- 对应 Surfe 文件/函数：`surfe_lib/lajaunie.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `_get_increment_pairs`、`process_input_data`、`get_interpolation_matrix`、`get_equality_values`、`get_inequality_values`、`setup_system_solver`、`convert_modified_kernel_to_rbf_kernel`、`eval_*`；`GRBF_Modelling_Methods::_update_interface_iso_values`。
- GeoRBF 目标模块：`src/model/lajaunie.rs`。
- 具体交付物：每层 reference point 和同层 difference 泛函；increment 索引/符号；普通、QP、restricted/modified/reconstruct 适用分支；iso-value 更新；标量/梯度评估。
- 必须运行的验证：多层、每层单/多点、重复/精确 level、inequality 和 restricted 数据的 C++ golden；reference/increment 离散快照；矩阵/RHS、残差、更新 level 和场；标准 fmt/clippy/test 全套。
- 完成门槛：reference point、pair 顺序和 iso-value 更新精确；全部源码可达纳入分支 parity；错误分类一致。
- 禁止事项：不得用全点绝对 value 代替 increment；不得容差合并 level；不得开始 Stratigraphic Horizons。

## T26 — Stratigraphic Horizons 层序与岩性约束

- 任务编号：T26。
- 目标：完成 Stratigraphic Horizons 的层序、岩性不等式、同层差值和最小层间约束。
- 依赖：T25。
- 对应 Surfe 文件/函数：`surfe_lib/stratigraphic_surfaces.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `_get_increment_pairs`、`_get_lithostratigraphic_increment_pairs_for_inequality_point`、`_get_closest_horizon_level_above_given_level`、`_get_closest_horizon_level_below_given_level`、`process_input_data`、`get_*matrix`、`get_*values`、`setup_system_solver`、`convert_modified_kernel_to_rbf_kernel`、`eval_*`。
- GeoRBF 目标模块：`src/model/stratigraphic.rs`。
- 具体交付物：精确层序；同层 difference；岩性点上下 horizon 配对；最小层间约束；普通/QP/restricted/reconstruct 可达路径；评估。
- 必须运行的验证：至少三层、边界层、层间/层外 inequality、相等和稀疏 level 的 C++ golden；pair、符号、layout、矩阵/RHS、可行性、层位关系和场；标准 fmt/clippy/test 全套。
- 完成门槛：层序、上下邻层、最小间隔、矩阵和最终层位关系 parity；错误与不可行分类一致。
- 禁止事项：不得排序成与源码相反的地层方向；不得臆造缺失层；不得开始 Continuous Property。

## T27 — Continuous Property 实际可达行为

- 任务编号：T27。
- 目标：只完成 Continuous Property 在冻结源码中实际可达的行为，并对部分实现、TODO 或不可达路径保持诚实分类。
- 依赖：T26。
- 对应 Surfe 文件/函数：`surfe_lib/continuous_property.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的构造/析构、`create_polynomial_basis`、`get_method_parameters`、`process_input_data`、`get_interpolation_matrix`、`get_equality_values`、`setup_system_solver`、`eval_*`、`measure_residuals`、`append_greedy_input`；对应 API 工厂分支。
- GeoRBF 目标模块：`src/model/continuous_property.rs`。
- 具体交付物：可达数据路径、矩阵、RHS、求解和评估；不可达/TODO 明确错误或文档；只对实际调用链存在的 Greedy 钩子预留内部归属，不宣称完成。
- 必须运行的验证：从公开 API 可达的 oracle 用例；矩阵/RHS、残差、标量/梯度；不可达/未实现分支的分类测试；标准 fmt/clippy/test 全套。
- 完成门槛：源码实际可达行为 parity；无能力臆增；TODO、空实现和不可达函数仍明确标识。
- 禁止事项：不得为对称性补写 C++ 不具备的模型能力；不得把类声明当完成；不得开始 Vector Field。

## T28 — Vector Field Planar Hessian、势函数与梯度

- 任务编号：T28。
- 目标：完成 Vector Field 的 Planar Hessian 系统、势函数标量评估和梯度输出。
- 依赖：T27。
- 对应 Surfe 文件/函数：`surfe_lib/vector_field.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的构造、`get_method_parameters`、`get_interpolation_matrix`、`get_equality_values`、`setup_system_solver`、`eval_scalar_interpolant_at_point`、`eval_vector_interpolant_at_point`；`RBFKernel::basis_planar_planar` 及 polynomial derivative 块。
- GeoRBF 目标模块：`src/model/vector_field.rs`。
- 具体交付物：Planar 分量 layout；Hessian 核系统；RHS 法向；polynomial derivative 约束；LU；势函数和向量梯度评估。
- 必须运行的验证：多方向 Planar 数据的 layout、Hessian matrix、RHS、残差、势函数和梯度 C++ golden；Hessian 对称/符号和有限差分势-梯度关系；标准 fmt/clippy/test 全套。
- 完成门槛：矩阵分量顺序、法向符号、势和梯度 parity；退化/缺少 Planar 的错误分类一致。
- 禁止事项：不得把梯度直接拟合作为不同问题；不得添加旋度/散度等新 API；不得开始最终公共 API。

### API、全局 parity 与性能

## T29 — 安全 Rust Builder 与 FittedModel

- 任务编号：T29。
- 目标：完成安全、类型化的 Rust Builder、不可变 `FittedModel`、单点和批量标量/梯度评估，同时复用已完成内部模型。
- 依赖：T28。
- 对应 Surfe 文件/函数：`surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的构造、全部 `Add*/Set*Constraints`、`Set*` 参数、`ComputeInterpolant`、`Get*Constraints`、`GetDataBoundsAndResolution`、`EvaluateInterpolantAtPoint(s)`、`EvaluateVectorInterpolantAtPoint(s)`、`GetNumberOfInterfaces`。
- GeoRBF 目标模块：`src/builder.rs`、`src/model/fitted.rs`、`src/lib.rs`。
- 具体交付物：可验证配置的 Builder；不可变、线程安全的 fitted model；五模型选择；单点/批量 scalar/gradient；输入所有权和错误；文档化生命周期；兼容层与惯用 Rust API 的边界。
- 必须运行的验证：五模型公共 API 端到端测试；单点/批量一致；拟合后不可变和并发只读测试；非法调用顺序和维数错误；文档测试；标准 fmt/clippy/test 全套。
- 完成门槛：所有已迁移可达功能能从安全公共 API 到达；批量与单点数值一致；无可变共享或 C++ 风格悬垂状态。
- 禁止事项：不得暴露内部矩阵可变引用；不得加入未迁移能力；不得通过 FFI 实现 API。

## T30 — Surfe 兼容名称、默认和错误分类测试

- 任务编号：T30。
- 目标：完成跨公共 API 的 Surfe 兼容名称、默认值、合法输入行为、拒绝行为和错误分类回归测试。
- 依赖：T29。
- 对应 Surfe 文件/函数：`modelling_parameters.h`、`grbf_exceptions.h`、`surfe_api.{h,cpp}`、`modeling_methods.cpp::get_method/create_rbf_kernel/check_input_data/check_interpolant@290dbe0ab344f4258a4935f05cad0f153f0f69a4`，以及五模型构造/检查分支。
- GeoRBF 目标模块：`tests/api_compat.rs`、`tests/error_compat.rs`、公共文档。
- 具体交付物：所有核/模型兼容名称表；逐字段默认快照；合法配置矩阵；拒绝/错误类别矩阵；调用顺序状态机测试；差异兼容文档。
- 必须运行的验证：全组合表驱动 oracle comparison；错误类别精确比较；默认构建端到端；标准 fmt/clippy/test/doc 全套。
- 完成门槛：离散 API 行为精确一致；安全修复有记录且不改变有效输入核心数学结果；没有未测试的公开配置分支。
- 禁止事项：不得用消息文本近似代替错误类别；不得扩大兼容名称；不得开始 Greedy 补全。

## T31 — Greedy 实际调用链审查与可达行为

- 任务编号：T31。
- 目标：审查 Greedy 的实际调用链，只实现冻结源码从纳入 API 实际可达的行为，不把 TODO、空实现或不可达代码宣称为完成。
- 依赖：T30。
- 对应 Surfe 文件/函数：`surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `run_greedy_algorithm`、`_output_greedy_debug_objects`；`modelling_input.{h,cpp}` 的 `Get_*_Indices_With_Large_Residuals`；各模型 `get_minimial_and_excluded_input`、`measure_residuals`、`append_greedy_input`；`Surfe_API::SetGreedyAlgorithm/ComputeInterpolant`。
- GeoRBF 目标模块：`src/greedy.rs`、各模型可达 Greedy 钩子、`docs/port/compatibility.md`。
- 具体交付物：调用链和可达性复核；实际可达的选择、残差测量、追加与停止行为；不可达/TODO 的明确文档和错误；确定性 trace。
- 必须运行的验证：对每个公开可达模型运行 oracle Greedy 数据；逐轮选中索引、残差、停止、最终矩阵和场 golden；不可达路径测试；标准 fmt/clippy/test 全套。
- 完成门槛：所有且仅有源码实际可达 Greedy 行为 parity；未实现/不可达部分未被误报；调试可视化输出保持排除。
- 禁止事项：不得补齐源码设想但未实现的 Greedy；不得迁移 debug visualization；不得创建新的优化任务。

## T32 — 全局差分 parity

- 任务编号：T32。
- 目标：执行全部模型、核、约束、求解器、错误路径和中间矩阵的全局差分 parity，并关闭所有既有兼容缺口。
- 依赖：T31。
- 对应 Surfe 文件/函数：T01 分类为迁移范围内且实际可达的全部 `surfe_lib/*` 与 `math_lib/*@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 符号；oracle 协议全部字段。
- GeoRBF 目标模块：`tests/parity/`、golden fixtures、`docs/port/parity-report.md`、必要的既有实现修正。
- 具体交付物：确定性正式 fixtures；全组合覆盖矩阵；离散 exact 比较；分层数值差分；三角验证汇总；缺陷兼容记录；按数据集的残差/可行性/预测报告。
- 必须运行的验证：完整 oracle fixture 再生成校验；全 Rust test suite；九核三角验证；五模型、普通/QP/LOQO/LU、errors、中间矩阵和 fields 差分；release 模式 smoke；fmt/clippy/doc/native 审计。
- 完成门槛：所有纳入范围的离散行为精确，数值行为达到既定分层容差；无未解释 mismatch、跳过或未运行测试；reference 身份精确；性能尚不宣称通过。
- 禁止事项：不得通过扩大容差、删 fixture、跳过失败或改变算法定义制造通过；不得做语义重构或性能近似；不得进入 T33 直到 parity 完整通过。

## T33 — Parity 保持下的性能优化

- 任务编号：T33。
- 目标：在 parity 不退化的前提下优化性能，并在同机、同数据、同线程、同优化条件下达到或超过冻结 Surfe。
- 依赖：T32。
- 对应 Surfe 文件/函数：各模型 `process_input_data/get_interpolation_matrix/setup_system_solver/eval_*`、`basis.cpp` 热路径、`matrix_solver.cpp` 与批量 API `Evaluate*AtPoints@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；以冻结构建作为基准而非迁移代码。
- GeoRBF 目标模块：`benches/`、`docs/port/performance-report.md`，及经 profile 证据支持的现有 Rust 热路径。
- 具体交付物：固定 benchmark 数据/查询；单线程和固定多线程 harness；预处理、装配、求解、标量评估、梯度评估、端到端分段报告；profile；不改变语义的优化及 parity 回归。
- 必须运行的验证：同机分别运行冻结 Surfe 和 GeoRBF release benchmark；同数据、查询点、线程数、优化级别并关闭调试/进度输出；每阶段多次测量并报告中位数；单线程和固定多线程；再次运行 T32 全 parity、fmt/clippy/test/native 审计。
- 完成门槛：所有 release-blocking 用例中 GeoRBF 中位时间不高于冻结 Surfe；不靠更多线程掩盖回退；六类时间都有结果；数值 parity 仍通过。无法测量或结果不确定即未完成并停留 T33。
- 禁止事项：未经实测不得声明性能达标；不得引入 FMM、低秩、局部近似、稀疏替代或新求解定义；不得牺牲 parity；不得只报告最好一次。

## T34 — 许可证、纯 Rust、平台与发布验收

- 任务编号：T34。
- 目标：完成许可证、来源追踪、纯 Rust 审计、平台 CI、文档和最终发布验收。
- 依赖：T33。
- 对应 Surfe 文件/函数：`License.txt`、`README.md`、全部实际迁移的 `path@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 来源清单，以及 CMake/Eigen/Qt/VTK/pybind 依赖边界；不新增算法对应项。
- GeoRBF 目标模块：根许可证/NOTICE（按审计结论）、`docs/`、CI、Cargo metadata/package、`docs/port/release-audit.md`。
- 具体交付物：逐模块来源追踪；许可证合规文件；纯 Rust 与发布包审计；Linux/macOS/Windows（按 Cargo 支持矩阵）CI 证据；用户 API 文档；parity 与性能报告链接；最终 release checklist。
- 必须运行的验证：`cargo fmt --check`、严格 clippy、全部 tests/doc tests、release build、`cargo package --list` 与打包试构建、`cargo tree --all-features`、native/build-script 二次审计、平台 CI、T32 parity 重跑、T33 benchmark 复核、`git diff --check`。
- 完成门槛：行为 parity、纯 Rust、性能不低于 Surfe、发布审计四项都有可复核证据并全部为真；发布包无 C++/reference/oracle；状态更新为 COMPLETE。
- 禁止事项：不得在任一证据缺失时完成；不得把 CI 未运行、reference 不可用或性能未判断写成通过；不得在 T34 后输出下一任务命令。

## 全局验收规则

以下规则约束 T00–T34，不得由单个任务放宽。

### 1. 行为等价

以下离散行为必须精确一致：

- 默认参数；
- 接受或拒绝的配置；
- 约束排序和去重；
- level 分组；
- reference point 和 increment pair；
- 约束自由度计数；
- 矩阵维数、行列顺序和分区；
- polynomial 项顺序；
- equality/inequality 的符号；
- 求解器分支；
- 成功或失败类别。

以下数值行为按 T04 固定的分层容差比较：

- 核值；
- 对第一个点和第二个点的一阶导数；
- 混合 Hessian；
- anisotropy transform；
- Modified Kernel；
- 完整矩阵和 RHS；
- LU/QP 残差与可行性；
- 标量场；
- 梯度场；
- QP 目标值和迭代证据。

权重向量不要求逐位相同；病态或非唯一系统应以残差、约束满足程度和预测场判断。不得仅因条件数很差而提前拒绝系统；必须按照冻结 Surfe 的实际求解语义尝试求解，再依据有限性、残差和约束可行性判断。

### 2. 三角验证

核函数和导数必须同时通过：

- 冻结 C++ Surfe golden comparison；
- 解析恒等式、符号和对称性检查；
- 有限差分数值检查。

Surfe 行为与数学检查冲突时不得静默修改。必须记录到既有兼容文档；有效输入的标准 GeoRBF 输出仍须满足既定 parity；API 安全修复不得改变核心数学结果。不复制 C++ 的未定义行为、内存问题或数据竞争；发现源码缺陷时，先在兼容测试中记录其可观察行为，再明确文档处理。

### 3. 纯 Rust

生产依赖中禁止出现：

- C++；
- Eigen；
- Qt；
- VTK；
- CMake runtime dependency；
- bindgen；
- `cxx`；
- 用于生产的 `cc` build dependency；
- MKL；
- OpenBLAS；
- LAPACK；
- native BLAS FFI。

外部 reference oracle 只可位于：

- `$SURFE_REFERENCE_DIR`；
- `../surfe`；
- `.cache/surfe-reference`；
- `.cache/surfe-oracle`。

reference 查找优先级为：环境变量 `SURFE_REFERENCE_DIR`，仓库同级目录 `../surfe`，被忽略的 `.cache/surfe-reference`。reference 必须校验 HEAD 或检出到固定提交。无法获得 reference 时，相关任务不得伪造通过。正常 Cargo 构建、测试和发布不得依赖 reference 或 oracle。

### 4. 范围与算法冻结

完成全局 parity 前，不进行会改变数学语义的算法重构，不引入 FMM、低秩、局部近似、稀疏替代或新的求解定义。现有 GeoRBF 代码必须先审查和复用；只有证据证明不兼容时才替换，不得无理由推倒重写。Qt、VTK、GUI、`geo_builder`、等值面显示和纯可视化代码始终排除。

不得把 TODO、空实现、不可达代码或类声明误报为 Surfe 已具备的能力。测试未执行、reference 无法运行、性能未判断或 parity 缺少证据时只能记录为未通过/未判断。

### 5. 性能

T33 只有满足以下条件才能完成：

- 同一机器；
- 同一数据和查询点；
- 同一线程数；
- 同一优化级别；
- 关闭调试和进度输出；
- 分别测量预处理、装配、求解、标量评估、梯度评估和端到端时间；
- 单线程和固定多线程均有结果；
- release-blocking 用例中，GeoRBF 中位运行时间不得高于冻结 Surfe；
- 不允许通过使用更多线程掩盖回退；
- 性能通过必须以数值 parity 已通过为前提；
- 无法判断只能记为未判断，不能记为通过。

若性能门槛未通过，T33 保持为当前任务，不能进入 T34。最终性能结论必须来自同机实测，未经实测不得声明达成。
