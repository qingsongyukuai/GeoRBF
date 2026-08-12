# T01 冻结 Surfe 字段级数据流

## 身份与边界

本数据流固定到 Surfe 提交 `290dbe0ab344f4258a4935f05cad0f153f0f69a4`。它描述有效输入如何从 `Surfe_API` 进入 `Constraints`，如何派生 layout、矩阵、RHS 和 weights，以及标量/梯度如何返回。本文不对 TODO、缺失定义或异常体做 T02 能力分类，也不把 Qt/VTK/geo_builder 数据流纳入迁移。

## 顶级字段流

```text
API 参数/矩阵
  -> Parameters + Constraints
  -> process_input_data 派生分组/pairs/bounds
  -> InternalParameters 固定计数与求解分支
  -> Kernel/Polynomial + constraint layout
  -> interpolation_matrix + equality/inequality RHS/范围
  -> System_Solver::weights
  -> eval_* 读取同一 layout 的 weights
  -> Point.scalar_field / Point._field_normal[3]
  -> API scalar/vector result
```

核心对象关系：

```text
Surfe_API
  owns/reference: method_: GRBF_Modelling_Methods*
  state: have_interpolant_, parameters_changed_, constraints_changed_

GRBF_Modelling_Methods / concrete model
  input: parameters, constraints
  derived: intern_params, interface_iso_values, interface_point_lists,
           interface_test_points, model-specific _increment_pairs
  basis: kernel, rbf_kernel, optional p_basis
  result: solver -> weights
```

冻结 C++ 使用裸指针并在部分复制/析构路径共享或别名这些对象。Rust 迁移应保持值和调用语义，但用所有权安全的值/枚举/共享只读表示，不能复制泄漏、悬垂或数据竞争。

## API 输入到约束字段

| API 输入 | 构造对象与字段 | 容器/标志 |
| --- | --- | --- |
| `AddInterfaceConstraint(x,y,z,level)` | `Point::_x/_y/_z/_c=0`；`Interface::_level=level,_residual=0,_level_bound={0,0}` | `constraints.itrface.push_back`；`parameters.use_interface=true`；`constraints_changed_=true` |
| `AddInequalityConstraint(x,y,z,level)` | `Inequality::_inequality_level=level,_residual=true` | `constraints.inequality`；`use_inequality=true` |
| `AddPlanarConstraintwNormal(x,y,z,nx,ny,nz)` | `_normal={nx,ny,nz}`，再派生 `_dip,_strike,_polarity` | `constraints.planar`；`use_planar=true` |
| `AddPlanarConstraintwStrikeDipPolarity` | `_dip,_strike,_polarity`，再派生 `_normal[3]` | `constraints.planar` |
| `AddPlanarConstraintwAzimuthDipPolarity` | `azimuth -> strike`，再走 strike/dip 构造 | `constraints.planar` |
| `AddTangentConstraint(x,y,z,tx,ty,tz)` | `_tangent={tx,ty,tz},_residual=0,_inner_product_constraint=0` | `constraints.tangent`；`use_tangent=true` |

四个 `Set*Constraints(MatrixXd)` 的列契约是 interface/inequality `x,y,z,level` 与 planar/tangent `x,y,z,vx,vy,vz`。每个 setter 先清空原容器，再调用 Add 路径；冻结 `SetTangentConstraints` 实际调用 planar Add 入口的事实保留到 T02/T30。

参数 setter 写入：

- `SetRBFKernel -> basis_type`；`SetRBFShapeParameter -> shape_parameter`；`SetPolynomialOrder -> polynomial_order`；`SetGlobalAnisotropy -> model_global_anisotropy`。
- `SetRestrictedRange -> use_restricted_range,interface_uncertainty,angular_uncertainty`。
- `SetRegressionSmoothing -> use_regression_smoothing,smoothing_amount`；冻结函数体对传入 bool 的具体赋值行为由 T06/T30 比较。
- `SetGreedyAlgorithm -> use_greedy,interface_uncertainty,angular_uncertainty`；具体 bool 行为同样由 T06/T30 比较。
- 所有 setter 令 `parameters_changed_=true`。

## 清洗、分组和派生输入

1. `remove_collocated_constraints` 对 `inequality -> itrface -> planar -> tangent` 各自以 `Point::operator<` 排序，再以 `collocated` 的三个坐标逐轴 `< Epilson` 判定相邻重复；类别之间不互相去重。
2. `get_interface_data` 清空并重建：
   - `interface_iso_values`：从 `Interface::level` 精确去重并降序；
   - `interface_point_lists[j]`：精确 level 相等的 Interface 副本；
   - `interface_test_points`：每组最先输入到分组的点；
   - 只有一个点的组仍先贡献 test point，随后从 `interface_point_lists` 删除。
3. Lajaunie `_increment_pairs`：每个保留组用 `[point_list[0], point_list[k+1]]`。
4. Stratigraphic `_increment_pairs` 按固定三段拼接：
   - 相邻 `interface_test_points[j], interface_test_points[j+1]`；
   - inequality 与最近的上/下 horizon test point；
   - 每层 reference 与同层其余点。
5. restricted range 派生：
   - Interface `level_bound={-interface_uncertainty,+interface_uncertainty}`；
   - Planar `normal_bound[axis][lower/upper]` 来自 strike/dip 角度界；
   - Tangent `angle_bound[lower/upper]` 来自角度 uncertainty。

这些容器都由具体模型持有，元素是原约束的值副本；更新 pair/test point 不自动更新所有原输入元素。

## 参数到核、多项式和求解分支

| 输入字段 | 派生/消费字段 | 消费节点 |
| --- | --- | --- |
| `model_type` | concrete `method_` | 两个工厂；公开 API 五分支，Greedy 工厂控制流见调用图 |
| `basis_type` | concrete `RBFKernel` | `create_rbf_kernel` |
| `shape_parameter` | Gaussian/MQ/MQ3/IMQ/WendlandC2/MaternC4 及 anisotropic 对应构造参数 | concrete kernel |
| `model_global_anisotropy` | ordinary/anisotropic factory branch | anisotropic 构造从 `constraints.planar` 派生 `_Transform` |
| `polynomial_order` | `n_poly_terms` 与 `Poly_Zero/First/Second` | model `get_method_parameters/create_polynomial_basis` |
| `use_inequality,use_restricted_range` 和约束实际计数 | `problem_type,modified_basis,poly_term,restricted_range` | model `get_method_parameters` |
| `use_regression_smoothing,smoothing_amount` | reference-point kernel value写入特定对角 | Single/Lajaunie matrix assembly |
| `min_stratigraphic_thickness` | sequenced horizon inequality RHS/range | Stratigraphic `get_inequality_values` |
| uncertainties | bounds、Greedy residual threshold | process/Greedy |

`setup_basis_functions` 总是保留 `rbf_kernel`；若 `intern_params.modified_basis`，`kernel` 指向新 `Modified_Kernel(rbf_kernel,interface_point_lists)`，否则 `kernel` 与 `rbf_kernel` 别名同一对象。

## 统一矩阵列含义

每个模型先由 `get_method_parameters` 写 `InternalParameters` 计数，再以相同顺序装配矩阵、RHS、weights 和评估。Planar 每点占三自由度 `{x,y,z}`，Tangent 每点占一自由度。

| 模型/阶段 | 核心列顺序 | 尺寸来源 |
| --- | --- | --- |
| Single Surface | `[inequality][interface][planar_x,planar_y,planar_z][tangent][poly]` | `n_constraints=n_ie+n_i+3n_p+n_t`；线性再加 `n_poly_terms` |
| Continuous Property 源码矩阵定义 | `[interface][planar_x,planar_y,planar_z][tangent][poly]` | 实际参数路径把 `n_p,n_t,n_poly` 置零；更宽定义只作 T02/T27 证据 |
| Vector Field | `[planar0_x,planar0_y,planar0_z,...]` | `3n_p`，无 poly |
| Lajaunie | `[same-level increment][planar_x,y,z][tangent][truncated poly]` | `_n_increment_pair+3n_p+n_t(+n_poly)` |
| Stratigraphic | `[sequenced interface increment][lithostrat inequality increment][same-level increment][planar_x,y,z][tangent][optional reconstructed poly]` | `_n_increment_pairs+3n_p+n_t(+n_poly)` |

### 核块字段流

- Value/Value -> `basis_pt_pt`。
- Value/Planar x/y/z -> `basis_pt_planar_x/y/z`；反向 -> `basis_planar_x/y/z_pt`。
- Planar/Planar -> 九个 `SecondDerivatives`，顺序 `DXDX,DXDY,DXDZ,DYDX,DYDY,DYDZ,DZDX,DZDY,DZDZ`。
- Value/Tangent、Tangent/Value -> `basis_pt_tangent,basis_tangent_pt`。
- Planar/Tangent 与反向 -> 三个 `FirstDerivatives`；Tangent/Tangent -> `basis_tangent_tangent`。
- increment `D(a,b)` 对 value center `c` 为 `K(a,c)-K(b,c)`；对另一个 increment `D(c,d)` 为 `K(a,c)-K(a,d)-K(b,c)+K(b,d)`。Planar/Tangent 块对 pair 两端使用相同相减规则。
- polynomial 下块 `P` 从约束位置的 `basis/dx/dy/dz` 得到；上块是转置 `P^T`；右下块为零。increment 模型使用两端 polynomial 值之差。

### RHS 字段流

- Single/Continuous interface equality读取 `Interface::level`；Planar 三项读取 `nx,ny,nz`；Tangent读取 `inner_product_constraint`（部分源码路径写常量零）；poly RHS 为零。
- Vector Field 每个 Planar 的 RHS 固定按 `nx,ny,nz` 排列。
- Lajaunie increment equality 为 `pair[0].level - pair[1].level`；Planar/Tangent/poly 后接。
- Stratigraphic ordinary QP：same-level pair进入 equality；相邻 horizon pair RHS 为 `min_stratigraphic_thickness`，lithostrat pair RHS 为零；Planar/Tangent 属 equality。
- Single ordinary inequality把每个 inequality level 映射为 `s(x)>=0` 或 `-s(x)>=0` 的行符号，RHS 零。
- restricted/LOQO 统一传 `b <= A*w <= b+r`：Interface/Planar/Tangent 的 lower/range 来自派生 bounds；Single lithology的范围还使用全点最大距离；Stratigraphic 三类 pair 分别使用 thickness、distance 和 uncertainty。

## 求解数据流

### Linear LU

```text
interpolation_matrix + equality_values
  -> Linear_LU_decomposition copies both
  -> Eigen partialPivLu().solve
  -> weights
  -> allFinite check
```

`validate_matrix_systems` 与 `check_solution` 是额外定义，不在模型直接拟合链上调用；T18/T02 保留其证据。

### 普通 predictor-corrector QP

```text
interpolation_matrix
  -> H = 2 * interpolation_matrix
equality rows/RHS -> A,b
inequality rows/RHS -> C,d
  -> Math_methods::quadratic_solver
       local x,y,z,s + KKT/residual/step state
  -> fvalues -> System_Solver::weights
```

### restricted/LOQO

```text
interpolation_matrix -> H = 2 * matrix
same matrix/selected A + lower b + range r
  -> Math_methods::quadratic_solver_loqo
       local x,y,g,z,t,s,v,w,p,q + predictor/corrector KKT
  -> w -> System_Solver::weights
```

求解对象构造时复制系统，模型随后只保留 `System_Solver* solver`。权重向量的索引必须与上述 layout 完全一致；后续 Rust 即使内部权重因病态/非唯一系统不同，也必须按全局规则比较残差、可行性和预测。

## weights 到标量/梯度字段

### Single/Continuous

对 evaluation `Point p`：

- scalar：按 layout 分段累加 `weight * basis_pt_pt/planar/tangent`，再加 `poly_basis * poly_weight`，写 `p.set_scalar_field`。
- gradient：value 列用 `basis_planar_{x,y,z}_pt`；Planar 列用 Hessian 三行；Tangent 列用 `basis_planar_tangent(DX/DY/DZ)`；poly 列用 `dx/dy/dz`，写 `p.set_vector_field`。

### Lajaunie/Stratigraphic

- scalar 的 pair 列为 `weight[k] * (K(p,pair0)-K(p,pair1))`，后接 Planar/Tangent/poly。
- gradient 对 pair 两端的一阶导数做同样相减，后接 Hessian/Tangent/poly。
- `_update_interface_iso_values` 以最终 weights 评估每个 `interface_test_points`，再用 `Point::scalar_field` 覆盖 `interface_iso_values`。

### Vector Field

- scalar/势函数：每个 Planar 的三 weights 乘 `basis_pt_planar_x/y/z`。
- gradient：每个 Planar 的三 weights 分别乘 Hessian 的 x/y/z 行，写 `Point::_field_normal[3]`。

API 单点评估还检查 API 状态位；批量路径把每行位置构造成独立 `Point` 并写入 Eigen 输出数组。批量标量使用 OpenMP，同时共享 `method_->kernel` 路径和进度计数；Rust 必须保留确定的数值语义而不能复制潜在数据竞争。

## QP 后重建数据流

Single/Lajaunie/Stratigraphic 的 `convert_modified_kernel_to_rbf_kernel` 使用当前 QP/Modified Kernel 场把不等式或 pair 转为确定的 equality 数据：

1. 在原约束或 pair 两端评估 scalar/gradient。
2. 把 scalar 写回新的/现有 Interface level，把 gradient 写回 Planar normal，把切向内积写回 Tangent。
3. `kernel = rbf_kernel`，关闭 modified/restricted，重算 equality 计数并启用 polynomial。
4. 重新装配普通 RBF 矩阵和 RHS，运行 Linear LU，替换 `solver`。
5. increment 模型再次更新 `interface_iso_values`。

这条路径不能在 Rust 中用 QP weights 直接跳过，因为最终评估读取的是重建后 layout 和 LU weights。

## 错误与状态数据流

- 低层 helper/模型以 `false` 或 `GRBF_Exceptions::*` 报错；模型 `setup_system_solver` 把 matrix/RHS/solver 阶段映射成不同异常常量。
- `ComputeInterpolant` 在 process/basis/solver 三段分别捕获 `std::exception`，用 `SurfeExceptions` 展平 nested exception 后重新抛出。
- `have_interpolant_` 在成功后为真；约束/参数修改令对应 changed 位为真；单点评估可抛 `missing_interpolant` 或 `interpolant_needs_update`。所有公共分支的精确行为由 T30 验证。
- 数组维数错误、未知核/模型、空间参数错误各有独立异常类型；Rust 不得只保留消息字符串。

## Eigen、CMake 与可视化边界

```text
CMake
  -> shared math_lib (Eigen)
  -> shared surfe_lib (math_lib + Eigen + optional OpenMP)
  -> surfepy (pybind11 + surfe_lib + math_lib)
  -> [GEO_BUILDER option only] Qt + VTK + geo_builder + test/main.cpp
```

- 迁移数据：参数、约束值、分组/pairs、核/多项式、矩阵/RHS、solver trace/weights、scalar/gradient、类型化错误。
- Rust 替代：Eigen matrix/vector、partial-pivot LU、LLT/特征分解、OpenMP 并行和 C++ exception/ownership。
- 排除：CMake 构建定义、DLL 宏、pybind ABI、Qt/VTK/GUI、geo_builder、regular grid/isosurface 写出与可视化、Windows console。
- 正常 Rust 生产/测试/发布不得依赖 Eigen、C++、CMake、pybind11、Qt、VTK、native BLAS/LAPACK/MKL/OpenBLAS、bindgen、`cxx` 或生产用 `cc`。

## 后续证据点

- T02 对本流中声明宽于实际调用链、空/常量/抛出体和状态差异做逐符号分类。
- T03 oracle 至少需要捕获清洗后约束、分组/pairs、layout 标签、kernel/derivative、matrix/RHS、solver 证据、weights、scalar/gradient 和 error category。
- T04 固定这些字段的离散精确比较与数值分层容差。
- T16–T21 固定 layout、装配、solver 与重建；T22–T28 分模型完成；T29/T30 固定安全 API；T31 只实现实际可达 Greedy。
