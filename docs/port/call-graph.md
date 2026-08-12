# T01 冻结 Surfe 调用图

## 身份、范围与读法

本图来自 `surfe_lib/*.{h,cpp}`、`math_lib/*.{h,cpp}`、`surfe_pybindings/pybindings.cpp` 和 `test/main.cpp` 在提交 `290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的声明、定义与调用点。箭头表示冻结源码中的直接调用或虚派发；“可能”只表示运行时参数分支，不表示 T02 已完成能力分类。

## 外部入口和对象工厂

```text
Python surfepy
  -> PYBIND11_MODULE(surfepy, m)
    -> Surfe_API 的同名构造/公共方法

C++ consumer
  -> Surfe_API(Parameters)
    -> Surfe_API::get_method_from_parameters
       Single_surface          -> new Single_Surface(params)
       Lajaunie_approach       -> new Lajaunie_Approach(params)
       Stratigraphic_horizons  -> new Stratigraphic_Surfaces(params)
       Vector_field            -> new Vector_Field(params)
       Continuous_property     -> new Continuous_Property(params)
       otherwise               -> unknown_modelling_mode

  -> Surfe_API(int)
       1 -> Single_Surface
       2 -> Lajaunie_Approach
       3 -> Vector_Field
       4 -> Stratigraphic_Surfaces
       5 -> Continuous_Property
       otherwise -> unknown_modelling_mode
```

`GRBF_Modelling_Methods::get_method(parameters)` 是 Greedy 内部的第二个工厂：Single/Lajaunie/Stratigraphic 分支显式构造对应类，其余值返回 `Continuous_Property`。它没有显式 Vector Field 分支。该控制流事实映射 T02/T31，不在 T01 判断其能力分类。

`test/main.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的活动路径是 `Geo_Builder::InitializeGRBFInterpolantObject(1) -> model.surfe->SetRestrictedRange -> LoadConstraintsFromFiles -> model.surfe->ComputeInterpolant`，然后进入 regular-grid、VTK 写出与可视化。前半段只证明核心入口被 geo_builder 调用；整个 harness 因 geo_builder/VTK/GUI 排除。

## 公共 API 状态与约束入口

```text
Add*/Set*Constraints
  -> 构造 Interface | Planar | Tangent | Inequality
  -> method_->constraints.{itrface|planar|tangent|inequality}
  -> parameters.use_* = true
  -> constraints_changed_ = true

参数 setter
  -> method_->parameters.<field>
  -> parameters_changed_ = true

Evaluate*AtPoint(s)
  -> 构造临时 Point
  -> method_->eval_scalar_interpolant_at_point
     或 method_->eval_vector_interpolant_at_point
  -> 读取 Point::{scalar_field|nx_interp,ny_interp,nz_interp}
```

矩阵 setter 先清空对应向量，再逐行调用 Add 入口。源码中 `SetTangentConstraints` 的直接调用目标是 `AddPlanarConstraintwNormal`；T01 保留该准确边，不替它改写意图。`GetDataBoundsAndResolution -> spatial_metrics`；`GetInterfaceReferencePoints` 读取 `method_->interface_test_points`；四类 getter 读取 `method_->constraints`。

## 顶级拟合调用链

```text
Surfe_API::ComputeInterpolant
  -> GRBF_Modelling_Methods::remove_collocated_constraints
       -> sort + unique(collocated)，分别处理四类约束
  -> virtual process_input_data
  -> virtual get_method_parameters
  -> GRBF_Modelling_Methods::setup_basis_functions
       -> create_rbf_kernel(parameters.basis_type, model_global_anisotropy)
       -> [modified_basis] new Modified_Kernel(rbf_kernel, interface_point_lists)
  -> virtual setup_system_solver
       -> virtual get_*matrix / get_*values
       -> System_Solver::solve
  -> have_interpolant_=true; constraints_changed_=parameters_changed_=false
```

`setup_basis_functions` 的普通核工厂分支：Cubic、Gaussian、IMQ、MQ、MQ3、R、TPS、WendlandC2、MaternC4；anisotropy 分支：ACubic、AGaussian、AIMQ、AMQ、AR、ATPS。六个 anisotropic 构造都调用 `RBFKernel::get_global_anisotropy(constraints.planar)`。

## 五模型从输入到求解

### Single Surface

```text
process_input_data
  -> get_interface_data
       -> _get_distinct_interface_iso_values
       -> _get_interface_points
  -> [restricted] Interface::setLevelBounds
  -> [restricted] Planar::setNormalBounds
  -> [restricted] Tangent::setAngleBounds

get_method_parameters
  -> layout [inequality | interface | planar x/y/z | tangent | optional poly]
  -> 普通 equality: Linear + ordinary kernel + polynomial
  -> inequality 或 restricted: Quadratic + Modified_Kernel + no polynomial

setup_system_solver
  Linear
    -> get_equality_values
    -> get_interpolation_matrix
    -> Linear_LU_decomposition::solve
  Quadratic, ordinary inequality
    -> get_inequality_values(vector), get_equality_values
    -> get_interpolation_matrix
    -> get_inequality_matrix + base get_equality_matrix
    -> Quadratic_Predictor_Corrector::solve
  Quadratic, restricted
    -> get_inequality_values(b,r)
    -> get_interpolation_matrix; inequality_matrix = interpolation_matrix
    -> Quadratic_Predictor_Corrector_LOQO::solve

convert_modified_kernel_to_rbf_kernel
  -> 用当前 QP 场更新 inequality/interface level、Planar normal、Tangent inner product
  -> kernel = rbf_kernel; inequality 清空；重建 layout/RHS/matrix
  -> Linear_LU_decomposition::solve
```

`get_interpolation_matrix` 的每个块调用 `Kernel::basis_pt_pt`、六个 point/planar 方向包装、`basis_planar_planar`、`basis_*_tangent`；普通路径再调用 `create_polynomial_basis -> _get_polynomial_matrix_block -> _insert_polynomial_matrix_blocks_in_interpolation_matrix`。smoothing 分支用两个参考点的 `kernel->basis_pt_pt()` 覆盖 inequality/interface 对角。

### Lajaunie Approach

```text
process_input_data
  -> get_interface_data
  -> _get_increment_pairs
       每个保留 level 组：reference point [0] 与其余点组成差值 pair
  -> [restricted] Planar/Tangent bounds

get_method_parameters
  -> layout [same-level increment | planar x/y/z | tangent | optional truncated poly]
  -> ordinary: Linear
  -> restricted: Quadratic Modified_Kernel/LOQO

setup_system_solver
  -> ordinary: equality/RBF matrix -> Linear_LU_decomposition::solve
  -> restricted: b/r + interpolation as A -> LOQO::solve
  -> _update_interface_iso_values

convert_modified_kernel_to_rbf_kernel
  -> 更新 pair 两端 level、Planar normal、Tangent inner product
  -> ordinary kernel + truncated poly + Linear LU
  -> _update_interface_iso_values
```

矩阵和评估中的 increment 泛函均展开为 `K(a,c)-K(a,d)-K(b,c)+K(b,d)`，point/planar/tangent 列分别做相同的两端差。`_update_interface_iso_values -> eval_scalar_interpolant_at_point(interface_test_points)`，再把结果写回 `interface_iso_values`。

### Stratigraphic Horizons

```text
process_input_data
  -> get_interface_data
  -> _get_increment_pairs
       1. 相邻 interface_test_points 的 sequenced interface pairs
       2. 每个 Inequality 与最近上/下 horizon 的 lithostratigraphic pairs
          -> _get_closest_horizon_level_above_given_level
          -> _get_closest_horizon_level_below_given_level
          -> Math_methods::sort_vector_w_index
       3. 各 level reference 与同层其余点的 equality pairs
  -> check_input_data
  -> [restricted] Planar/Tangent bounds

get_method_parameters
  -> layout [全部 increment | planar x/y/z | tangent]
  -> ordinary QP: 前两类 pair 是 inequality，第三类 pair + planar/tangent 是 equality
  -> restricted LOQO: 每个 layout 行为区间约束

setup_system_solver
  -> ordinary: get_inequality_values + get_equality_values
       -> interpolation/equality/inequality matrices
       -> Quadratic_Predictor_Corrector::solve
  -> restricted: get_inequality_values(b,r)
       -> interpolation as A -> LOQO::solve
  -> _update_interface_iso_values

convert_modified_kernel_to_rbf_kernel
  -> 更新所有 pair/Planar/Tangent
  -> 把全部 pair 改作 equality，ordinary kernel + truncated poly
  -> Linear_LU_decomposition::solve
  -> _update_interface_iso_values
```

### Continuous Property

```text
process_input_data
  -> 要求 constraints.itrface 非空
get_method_parameters
  -> n_interface = itrface.size
  -> n_inequality=n_planar=n_tangent=0; n_poly_terms=0
  -> Linear, ordinary kernel
setup_system_solver
  -> get_equality_values(level)
  -> get_interpolation_matrix
  -> Linear_LU_decomposition::solve
eval_scalar/vector
  -> 通用 interface/planar/tangent 核和可选 poly 读取路径
```

类中仍声明或定义更宽的 matrix、residual、append 与转换槽；T01 把所有调用边和定义列入清单，不把声明宽度误报成公共可达能力。其最终归类在 T02/T27/T31。

### Vector Field

```text
process_input_data -> 内联空体（事实记录，分类在 T02）
get_method_parameters
  -> n_planar = constraints.planar.size; n_constraints=n_equality=3*n_planar
  -> no polynomial, ordinary kernel, Linear
setup_system_solver
  -> get_equality_values [nx,ny,nz] per Planar
  -> get_interpolation_matrix [3P x 3P Hessian blocks]
  -> Linear_LU_decomposition::solve
eval_scalar
  -> Σ weights[3k..3k+2] * basis_pt_planar_{x,y,z}
eval_vector
  -> Σ weights * basis_planar_planar(九个 Hessian 分量)
```

## 核内部调用图

```text
model matrix/evaluation
  -> Kernel::set_points(p1,p2)
  -> Kernel::basis_* virtual combination

ordinary RBFKernel combination
  basis_pt_pt          -> basis
  point/planar         -> dx/dy/dz on p1 or p2
  planar/planar        -> one of dxx..dzz
  point/tangent        -> p2 tangent dot [dx_p2,dy_p2,dz_p2]
  tangent/point        -> p1 tangent dot [dx_p1,dy_p1,dz_p1]
  tangent/tangent      -> t1^T Hessian t2
  planar/tangent       -> selected Hessian row dot t2
  tangent/planar       -> selected Hessian column dot t1

isotropic concrete kernel
  -> radius -> Point::{x,y,z,c}
anisotropic concrete kernel
  -> scaled_radius -> Matrix3f _Transform * spatial delta

Modified_Kernel combination
  -> Lagrangian_Polynomial_Basis::{poly,poly_dx,poly_dy,poly_dz}
  -> _aRBFKernel 在 evaluation/unisolvent/center 组合上求值或导数
  -> 消去式 base - t1 - t2 + t3 + t4
```

`Lagrangian_Polynomial_Basis` 构造调用 `_get_unisolvent_subset -> Math_methods::sort_vector_w_index -> _initialize_basis`。`Poly_Zero/First/Second` 由模型的 `create_polynomial_basis` 创建；模型先 `set_point`，再调用 `basis/dx/dy/dz`。

## 求解器内部调用图

```text
Linear_LU_decomposition::solve
  -> Eigen::partialPivLu().solve
  -> weights.allFinite

Quadratic_Predictor_Corrector::solve
  -> Math_methods::quadratic_solver(H=2*interpolation,Aeq,Cineq,b,d)
       -> 多次 Eigen::partialPivLu().solve(KKT,rhs)
       -> _find_step_length
       -> max_element_wrt_zero

Quadratic_Predictor_Corrector_LOQO::solve
  -> Math_methods::quadratic_solver_loqo(H,A,b,r)
       -> 初始/预测/校正 KKT partialPivLu
       -> _find_positivity_step -> _find_step
```

`validate_matrix_systems` 和 `Linear_LU_decomposition::check_solution` 有定义，但模型 `setup_system_solver` 的直接路径只调用各 solver 的 `solve`。是否可达、部分实现或缺陷由 T02 分类。

## Greedy 调用链

```text
GRBF_Modelling_Methods::run_greedy_algorithm
  -> 参数 uncertainty 非零检查
  -> get_method(parameters) 生成 greedy_method
  -> greedy_method->constraints.compute_avg_nn_distances
       -> 四类 compute_* -> avg_nn_distance
  -> this->get_minimial_and_excluded_input
  -> greedy_method->constraints = greedy_input
  -> loop
       greedy_method->process_input_data
       greedy_method->get_method_parameters
       greedy_method->setup_basis_functions
       greedy_method->setup_system_solver
       greedy_method->measure_residuals(excluded_input)
       greedy_method->append_greedy_input(excluded_input)
       greedy_method->_SetIteration
```

Single/Lajaunie 的候选追加调用四个 `Get_*_STL_Vector_Indices_With_Large_Residuals`，进而调用 `distance_btw_pts` 与 `Math_methods::sort_vector_w_index`。Continuous Property 有自己的优先追加定义。Stratigraphic/Vector 的相应槽体、`_output_greedy_debug_objects` 与注释中的调试输出只记录为源码边界；T31 只能实现 T02 证明的实际可达行为。

## 调用者、被调用者和状态所有者核对

| 节点 | 主要调用者 | 主要被调用者 | 状态所有者 |
| --- | --- | --- | --- |
| `Surfe_API` | C++/pybind consumer | 模型虚接口、空间函数 | `Surfe_API` 拥有 API 状态位和 `method_` |
| `GRBF_Modelling_Methods` | `Surfe_API`、Greedy | 工厂、核、具体模型虚派发 | 具体模型对象拥有约束、参数、派生分组和指针 |
| 具体模型 `process/get/setup/eval` | API/Greedy | Kernel、Polynomial、Solver、共享 helper | 具体模型拥有 layout 计数与 pair 状态 |
| `Kernel`/`RBFKernel`/`Modified_Kernel` | 模型装配和评估 | Point、Lagrangian、具体导数 | kernel 对象暂存 `_p1/_p2` 和核参数 |
| `System_Solver` | 模型 `setup_system_solver` | Math_methods/Eigen | solver 复制输入系统并拥有 `weights` |
| `Math_methods` | solver、空间/Greedy/model pair helper | Eigen/STL | 无持久对象状态；LOQO/PC 状态为函数局部 |

## 后续归属

- T02 关闭所有定义/声明、TODO、抛出体、常量体和可达性分类问题。
- T03/T04 依据本图选择 API、矩阵、RHS、weights、fields、errors 和 solver trace oracle 点。
- T16/T17 固定五模型 layout 与装配；T18–T21 固定求解/重建；T22–T28 分模型实现。
- T29/T30 固定安全 API 和状态机；T31 复核并实现实际可达 Greedy。未创建新的顶级任务。
