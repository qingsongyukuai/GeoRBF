# T01 冻结 Surfe C++ 清单

## 审计身份与记法

- 权威源码：`https://github.com/MichaelHillier/surfe.git`，提交 `290dbe0ab344f4258a4935f05cad0f153f0f69a4`。
- 本清单覆盖 `surfe_lib/*.{h,cpp}`、`math_lib/*.{h,cpp}`、`surfe_pybindings/pybindings.cpp` 到核心的适配，以及 `test/main.cpp`。
- `geo_builder`、Qt、VTK、GUI、等值面和纯可视化保持排除；它们只在依赖边界与 `test/main.cpp` 的调用证据中出现。
- 下文 `Type::{a,b}` 是逐符号记法，明确展开为 `Type::a` 与 `Type::b` 两个独立符号；同理，公共方法集 `K` 用于避免十五次复制相同的 17 个虚函数名。
- 本任务只清点声明、定义和调用证据。TODO、常量返回、抛出、缺失定义、可达性与缺陷的唯一能力分类属于 T02，不能从本清单推断“已实现”。

## 逐文件清单

| 冻结文件 | 范围内内容 | 依赖处置 |
| --- | --- | --- |
| `surfe_lib/surfe_api.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | `Surfe_API` 公共构造、约束、参数、拟合、标量/向量评估和查询入口 | 迁移，T06/T29/T30 |
| `surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 五模型工厂、API 状态机、批量 OpenMP 评估 | 迁移；OpenMP 用纯 Rust 并发替代，T29/T30 |
| `surfe_lib/modelling_parameters.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 常量、枚举、`Parameters`/`InternalParameters`/`InputParameters` | 迁移，T06 |
| `surfe_lib/modelling_input.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 约束值类型、内联访问器、空间函数声明 | 迁移，T07–T09/T31 |
| `surfe_lib/modelling_input.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 方向转换、空间统计、Greedy 残差候选筛选 | 迁移，T07/T09/T31 |
| `surfe_lib/modeling_methods.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 抽象模型基类、共享状态、全部虚函数槽 | 迁移，T08/T16/T29/T31 |
| `surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 分组/去重、核和模型工厂、共享拟合与 Greedy 流程 | 迁移，T08/T12–T14/T31 |
| `surfe_lib/basis.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 核、各向异性核、多项式、Lagrangian 与 Modified Kernel 类型树 | 迁移，T10–T15 |
| `surfe_lib/basis.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 半径/变换、核与导数、核组合、多项式和 Modified Kernel 定义 | 迁移，T10–T15；Eigen 特征分解用纯 Rust 替代 |
| `surfe_lib/matrix_solver.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 三个求解器对象及其矩阵/向量所有权 | 迁移，T18–T20 |
| `surfe_lib/matrix_solver.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | LU、普通 QP 与 LOQO 调用桥 | 迁移，T18–T20；Eigen 求解用纯 Rust 替代 |
| `surfe_lib/single_surface.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | `Single_Surface` 覆写和私有装配辅助 | 迁移，T22–T24 |
| `surfe_lib/single_surface.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | Single Surface 全路径 | 迁移，T22–T24 |
| `surfe_lib/lajaunie.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | `Lajaunie_Approach` 覆写、increment 状态 | 迁移，T25 |
| `surfe_lib/lajaunie.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | reference point、同层 increment、求解与评估 | 迁移，T25 |
| `surfe_lib/stratigraphic_surfaces.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | `Stratigraphic_Surfaces` 覆写、层序计数 | 迁移，T26 |
| `surfe_lib/stratigraphic_surfaces.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 层序/岩性 increment、QP、重建和评估 | 迁移，T26 |
| `surfe_lib/continuous_property.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | `Continuous_Property` 覆写 | 证据限定迁移，T02/T27 |
| `surfe_lib/continuous_property.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | property equality、矩阵、求解、评估与 Greedy 钩子 | 证据限定迁移，T02/T27/T31 |
| `surfe_lib/vector_field.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | `Vector_Field` 覆写 | 证据限定迁移，T02/T28 |
| `surfe_lib/vector_field.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | Planar Hessian、LU、势与梯度 | 迁移，T28 |
| `surfe_lib/grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 22 个异常类别、嵌套异常扁平化和常量实例 | 迁移，T06/T30 |
| `surfe_lib/debug.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | Windows console 声明；其余模板为注释 | 排除可视化/调试，T02 确认 |
| `surfe_lib/debug.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | `open_console_window` Windows 实现 | 排除 |
| `surfe_lib/surfe_lib_module.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | DLL export/deprecation 宏 | C++ ABI 构建边界，Rust 不迁移 |
| `math_lib/math_methods.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | QP、角度、排序和步长辅助声明 | 迁移，T08/T19/T20/T31 |
| `math_lib/math_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 普通/LOQO QP 算法、排序与角度定义 | 迁移，T08/T19/T20/T31 |
| `math_lib/math_lib_module.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | DLL export/deprecation 宏 | C++ ABI 构建边界，Rust 不迁移 |
| `surfe_pybindings/pybindings.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | `PYBIND11_MODULE(surfepy, m)` 与 30 个 `Surfe_API` 适配绑定 | 仅作为公开可达性证据；pybind11 排除 |
| `test/main.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | `main` 经 `Geo_Builder::surfe` 调用 restricted range 与 `ComputeInterpolant` | harness 依赖 geo_builder/VTK/GUI，排除；核心调用只作证据 |

## 类型、常量与字段符号

### 参数和输入类型

- 宏：`D2R`、`R2D`、`Epilson`。
- `Parameter_Types::{DWRT,SecondDerivatives,FirstDerivatives,RBF,SolverType,ModelType,AXIS}`；枚举值依次为 `{PT1,PT2}`、`{DXDX,DXDY,DXDZ,DYDX,DYDY,DYDZ,DZDX,DZDY,DZDZ}`、`{DX,DY,DZ}`、`{Cubic,Gaussian,MQ,MQ3,IMQ,TPS,R,WendlandC2,MaternC4}`、`{Linear,Quadratic}`、`{Single_surface,Lajaunie_approach,Stratigraphic_horizons,Continuous_property,Vector_field}`、`{Xaxis,Yaxis,Zaxis}`。
- `Parameters` 的字段：`model_type`、`min_stratigraphic_thickness`、`use_interface`、`use_planar`、`use_tangent`、`use_inequality`、`basis_type`、`shape_parameter`、`polynomial_order`、`advanced_parameters`、`model_global_anisotropy`、`use_greedy`、`use_restricted_range`、`smoothing_amount`、`use_regression_smoothing`、`interface_uncertainty`、`angular_uncertainty`，以及 `Parameters::Parameters`。
- `InternalParameters` 的字段：`n_interface`、`n_planar`、`n_inequality`、`n_tangent`、`n_constraints`、`n_equality`、`modified_basis`、`poly_term`、`n_poly_terms`、`problem_type`、`restricted_range`，以及 `InternalParameters::InternalParameters`。
- `InputParameters::{parameters,interface_file,planar_file,tangent_file,inequality_file}`。

### 约束值类型

- `Point` 字段：`_x,_y,_z,_c,_scalar_field,_field_normal[3]`；符号 `Point::{Point,set_x,set_y,set_z,set_c,scalar_field,set_scalar_field,set_vector_field,operator<,x,y,z,c,nx_interp,ny_interp,nz_interp}` 及自由函数 `collocated`。
- `Interface` 字段：`_level,_residual,_level_bound[2]`；符号 `Interface::{Interface,level,residual,level_lower_bound,level_upper_bound,setResidual,setLevel,setLevelBounds}`。
- `Inequality` 字段：`_inequality_level,_residual`；符号 `Inequality::{Inequality,level,residual,setResidual}`。
- `Planar` 字段：`_dip,_strike,_polarity,_normal[3],_residual,_normal_bound[3][2]`；符号 `Planar::{Planar,_compute_strike_dip_polarity_from_normal,_compute_normal_from_strike_dip_polarity,getDipVector,getStrikeVector,dip,strike,polarity,nx,ny,nz,nx_lower_bound,nx_upper_bound,ny_lower_bound,ny_upper_bound,nz_lower_bound,nz_upper_bound,setNormalBounds,residual,setResidual,setNormal}`。
- `Tangent` 字段：`_tangent[3],_residual,_angle_bound[2],_inner_product_constraint`；符号 `Tangent::{Tangent,tx,ty,tz,residual,angle_lower_bound,angle_upper_bound,inner_product_constraint,setResidual,setAngleBounds,setInnerProductConstraint}`。
- `Constraints` 拥有四个向量 `inequality,itrface,planar,tangent` 和四个平均最近邻字段；符号 `Constraints::{Constraints,compute_inequality_avg_nn_distance,compute_interface_avg_nn_distance,compute_planar_avg_nn_distance,compute_tangent_avg_nn_distance,compute_avg_nn_distances,GetInequalityAvgNNDist,GetInterfaceAvgNNDist,GetPlanarAvgNNDist,GetTangentAvgNNDist,SetInequalityAvgNNDist,SetInterfaceAvgNNDist,SetPlanarAvgNNDist,SetTangentAvgNNDist}`。
- `SpatialParameters::{resolution,xmin,xmax,ymin,ymax,zmin,zmax}`。
- 自由空间符号：`convert_constraints_to_points`、`distance_btw_pts`、`nearest_neighbour_index`、`get_n_nearest_neighbours_to_point`、`furtherest_neighbour_index(Point,points)`、`furtherest_neighbour_index(points,points)`、`avg_nn_distance`、`spatial_metrics`、`Find_STL_Vector_Indices_FurtherestTwoPoints`、`Find_STL_Vector_Index_ofPointClosestToOtherPointWithinDistance`、`calculate_bounds`、`get_extremal_point_data_indices_from_points`、`is_index_in_list`、`get_largest_distance_between_points`、`Get_Inequality_STL_Vector_Indices_With_Large_Residuals`、`Get_Interface_STL_Vector_Indices_With_Large_Residuals`、`Get_Planar_STL_Vector_Indices_With_Large_Residuals`、`Get_Tangent_STL_Vector_Indices_With_Large_Residuals`、`get_maximal_axial_variability_order`。

## 核与多项式逐符号清单

### 抽象层和组合层

- `Kernel` 字段 `_p1,_p2`；符号 `Kernel::{~Kernel,set_points,set_center,set_evaluation_points,p1,p2,basis_pt_pt,basis_pt_planar_x,basis_planar_x_pt,basis_pt_planar_y,basis_planar_y_pt,basis_pt_planar_z,basis_planar_z_pt,basis_pt_tangent,basis_tangent_pt,basis_planar_planar,basis_tangent_tangent,basis_planar_tangent,basis_tangent_planar,clone}`。
- `RBFKernel` 字段 `_radius,_x_delta,_y_delta,_z_delta,_c_delta,_Global_Plunge[3],_Transform`；符号 `RBFKernel::{RBFKernel,~RBFKernel,radius,scaled_radius,get_global_anisotropy,basis,dx_p1,dx_p2,dy_p1,dy_p2,dz_p1,dz_p2,dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz,basis_pt_pt,basis_pt_planar_x,basis_planar_x_pt,basis_pt_planar_y,basis_planar_y_pt,basis_pt_planar_z,basis_planar_z_pt,basis_pt_tangent,basis_tangent_pt,basis_planar_planar,basis_tangent_tangent,basis_planar_tangent,basis_tangent_planar,clone}`。
- 每个具体核的公共方法集 `K` 精确为 `{basis,dx_p1,dx_p2,dy_p1,dy_p2,dz_p1,dz_p2,dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz,clone}`。
- `Cubic::{~Cubic,K}`、`Gaussian::{Gaussian,~Gaussian,K}`、`MQ::{MQ,~MQ,K}`、`MQ3::{MQ3,~MQ3,K}`、`TPS::{~TPS,K}`、`IMQ::{IMQ,~IMQ,K}`、`R::{~R,K}`、`WendlandC2::{WendlandC2,~WendlandC2,K}`、`MaternC4::{MaternC4,~MaternC4,K}`。
- `ACubic::{ACubic,~ACubic,K}`、`AGaussian::{AGaussian,~AGaussian,K}`、`AMQ::{AMQ,~AMQ,K}`、`ATPS::{ATPS,~ATPS,K}`、`AIMQ::{AIMQ,~AIMQ,K}`、`AR::{AR,~AR,K}`。这些六类经 `RBFKernel::get_global_anisotropy` 填充 `Matrix3f _Transform`。
- `Modified_Kernel` 拥有 `_aRBFKernel` 与 `_aLPB`；符号 `Modified_Kernel::{Modified_Kernel,Modified_Kernel(copy),~Modified_Kernel,basis_pt_pt,basis_pt_planar_x,basis_planar_x_pt,basis_pt_planar_y,basis_planar_y_pt,basis_pt_planar_z,basis_planar_z_pt,basis_pt_tangent,basis_tangent_pt,basis_planar_planar,basis_tangent_tangent,basis_planar_tangent,basis_tangent_planar,clone}`。

### 多项式

- `Polynomial_Basis` 字段 `_p,_truncated`；符号 `Polynomial_Basis::{Polynomial_Basis,set_point,basis,dx,dy,dz,clone}`。
- `Poly_Zero::{Poly_Zero,basis,dx,dy,dz,clone}`、`Poly_First::{Poly_First,basis,dx,dy,dz,clone}`、`Poly_Second::{Poly_Second,basis,dx,dy,dz,clone}`。
- `Lagrangian_Polynomial_Basis` 字段 `_polynomial_constants,_derivative_polynomial_constants,unisolvent_subset_points`；符号 `Lagrangian_Polynomial_Basis::{Lagrangian_Polynomial_Basis,_get_unisolvent_subset,_initialize_basis,poly,poly_dx,poly_dy,poly_dz}`。

## 模型基类、五模型和覆写清单

### 共享基类

- `GRBF_Modelling_Methods` 拥有 `intern_params,_iteration,interface_iso_values,interface_point_lists,constraints,parameters,solver,kernel,rbf_kernel,error_msg,interface_test_points`。
- 全部符号：`GRBF_Modelling_Methods::{~GRBF_Modelling_Methods,_get_distinct_interface_iso_values,_get_interface_points,_get_distinct_inequality_iso_values,_interface_points_are_coplanar,get_interface_data,check_input_data,_update_interface_iso_values,_output_greedy_debug_objects,_SetIteration,get_method,create_rbf_kernel,get_interface_points_ouput,remove_collocated_constraints,get_interface_iso_values,setup_basis_functions,check_interpolant,run_greedy_algorithm,get_equality_matrix,get_interpolation_matrix,get_equality_values,eval_scalar_interpolant_at_point,eval_vector_interpolant_at_point,get_method_parameters,process_input_data,setup_system_solver,get_minimial_and_excluded_input,measure_residuals,append_greedy_input,convert_modified_kernel_to_rbf_kernel,clone}`。
- 十二个必须由具体模型归属的槽位 `V`：`{get_interpolation_matrix,get_equality_values,eval_scalar_interpolant_at_point,eval_vector_interpolant_at_point,get_method_parameters,process_input_data,setup_system_solver,get_minimial_and_excluded_input,measure_residuals,append_greedy_input,convert_modified_kernel_to_rbf_kernel,clone}`。

### 覆写矩阵

| 模型 | 构造/私有状态 | `V` 覆写归属与额外符号 |
| --- | --- | --- |
| `Single_Surface` | `Single_Surface::{Single_Surface,~Single_Surface,_get_polynomial_matrix_block,_insert_polynomial_matrix_blocks_in_interpolation_matrix}`，拥有 `p_basis` | `Single_Surface::V`；额外 `create_polynomial_basis,get_inequality_matrix,get_inequality_values(vector),get_inequality_values(b,r)` |
| `Lajaunie_Approach` | `Lajaunie_Approach::{Lajaunie_Approach,~Lajaunie_Approach,_get_increment_pairs,_get_polynomial_matrix_block,_insert_polynomial_matrix_blocks_in_interpolation_matrix}`，拥有 `_n_increment_pair,_increment_pairs,p_basis` | `Lajaunie_Approach::V`；额外 `create_polynomial_basis,get_inequality_values(b,r)` |
| `Stratigraphic_Surfaces` | `Stratigraphic_Surfaces::{Stratigraphic_Surfaces,~Stratigraphic_Surfaces,_get_increment_pairs,_get_lithostratigraphic_increment_pairs_for_inequality_point,_get_closest_horizon_level_above_given_level,_get_closest_horizon_level_below_given_level,_get_polynomial_matrix_block,_insert_polynomial_matrix_blocks_in_interpolation_matrix}`，拥有四个 pair 计数、`_increment_pairs,p_basis` | `Stratigraphic_Surfaces::V`；额外 `create_polynomial_basis,get_inequality_matrix,get_inequality_values(vector),get_inequality_values(b,r)`；三个 Greedy 槽的内联体记录于源码，能力分类留给 T02 |
| `Continuous_Property` | `Continuous_Property::{Continuous_Property,~Continuous_Property,_get_polynomial_matrix_block,_insert_polynomial_matrix_blocks_in_interpolation_matrix}`，拥有 `p_basis` | `Continuous_Property::V`；额外 `create_polynomial_basis`；`get_minimial_and_excluded_input` 与 `convert_modified_kernel_to_rbf_kernel` 的内联体记录于源码，分类留给 T02 |
| `Vector_Field` | `Vector_Field::{Vector_Field,~Vector_Field}` | `Vector_Field::V`；`process_input_data`、三个 Greedy 槽和转换槽的内联体记录于源码，分类留给 T02 |

工厂审计：`Surfe_API::get_method_from_parameters` 明列五个 `ModelType` 分支；`Surfe_API(int)` 明列整数 1–5；`GRBF_Modelling_Methods::get_method` 明列 Single/Lajaunie/Stratigraphic 并以 `Continuous_Property` 为默认返回，没有显式 `Vector_Field` 分支。这里仅记录控制流差异，T02/T31 决定能力与可达性分类。

## 求解、数学、API 和异常逐符号清单

### 求解器与数学辅助

- `System_Solver`：`weights`；`System_Solver::{System_Solver,~System_Solver,solve,validate_matrix_systems}`。
- `Linear_LU_decomposition`：`_interpolation_matrix,_constraint_values`；`Linear_LU_decomposition::{Linear_LU_decomposition,~Linear_LU_decomposition,solve,validate_matrix_systems,check_solution}`。
- `Quadratic_Predictor_Corrector`：`_interpolation_matrixD,_hessian_matrixD,_equality_matrixD,_inequality_matrixD,_equality_vectorD,_inequality_vectorD`；`Quadratic_Predictor_Corrector::{Quadratic_Predictor_Corrector,solve,validate_matrix_systems}`。
- `Quadratic_Predictor_Corrector_LOQO`：`_H,_A,_b,_r`；`Quadratic_Predictor_Corrector_LOQO::{Quadratic_Predictor_Corrector_LOQO,solve,validate_matrix_systems}`。
- `Math_methods::{_find_step_length,_rot,_get_double,_find_step,_find_positivity_step,sort_vector_w_index,max_element_wrt_zero,SWAP,angle_btw_2_vectors,RandomDouble,quadratic_solver,quadratic_solver_loqo}`。声明/定义双向差异（包括 `_rot`）留给 T02 唯一分类。

### `Surfe_API`

- 私有字段：`method_,have_interpolant_,parameters_changed_,constraints_changed_`。
- 符号：`Surfe_API::{progress,get_method_from_parameters,Surfe_API,AddInterfaceConstraint,AddPlanarConstraintwNormal,AddPlanarConstraintwStrikeDipPolarity,AddPlanarConstraintwAzimuthDipPolarity,AddTangentConstraint,AddInequalityConstraint,ComputeInterpolant,SetRegressionSmoothing,SetGreedyAlgorithm,SetRestrictedRange,SetRBFKernel(enum),SetRBFKernel(name),SetRBFShapeParameter,SetPolynomialOrder,SetGlobalAnisotropy,EvaluateInterpolantAtPoint,EvaluateInterpolantAtPoints,EvaluateVectorInterpolantAtPoint,EvaluateVectorInterpolantAtPoints,GetDataBoundsAndResolution,GetInterfaceReferencePoints,GetInterfaceConstraints,SetInterfaceConstraints,GetPlanarConstraints,SetPlanarConstraints,GetTangentConstraints,SetTangentConstraints,GetInequalityConstraints,SetInequalityConstraints,GetNumberOfInterfaces,InterpolantComputed}`。
- `surfe_pybindings/pybindings.cpp` 绑定除 enum 重载、`GetNumberOfInterfaces` 与 `InterpolantComputed` 外的上述公共构造/方法；这是 Python 到相同 `Surfe_API` 对象的适配，不拥有算法状态。

### 异常

- 异常类及各自 `what`：`nointerfacedata`、`nointerfaceincrementpairs`、`noplanardata`、`invalidinputdata`、`failurecomputingglobalanisotropy`、`failurecreatinganisotropickernel`、`failuresettingupbasisfunctions`、`failurecreatingmodifiedkernel`、`failurecreatinglagrangianpolynomialbasis`、`linearsolverfailure`、`pcquadratricsolverfailure`、`loqoquadratricsolverfailure`、`errorcomputinginterpolationmatrix`、`errorcomputingequalityvector`、`errorcomputinginequalityvector`、`errorupdatinginterfaceisovalues`、`errorcomputinginterpolant`、`missinginterpolant`、`unknownrbf`、`interpolantneedsupdate`、`unknownmodellingmode`、`problemcomputingspatialparameters`、`arrayhasincorrectdimensions`。
- `SurfeExceptions::{append_exceptions,SurfeExceptions,what}`。
- `GRBF_Exceptions` 常量实例逐一为 `no_iterface_data,no_interface_increment_pairs,no_planar_data,invalid_input_data,failure_computing_global_anisotropy,failure_creating_anisotropic_kernel,failure_setting_up_basis_functions,failure_creating_modified_kernel,failure_creating_lagrangian_polynomial_basis,linear_solver_failure,pc_quadratic_solver_failure,loqo_quadratic_solver_failure,error_computing_interpolation_matrix,error_computing_equality_vector,error_computing_inequality_vector,error_updating_interface_iso_values,error_computing_interpolant,missing_interpolant,unknown_rbf,interpolant_needs_update,unknown_modelling_mode,problem_computing_spatial_parameters,array_has_incorrect_dimensions`。

### 构建/适配符号

- `open_console_window` 仅在 `_WIN32` 分支分配 console，排除。
- `SURFE_LIB_EXPORT`/`SURFE_LIB_NO_EXPORT`/`SURFE_LIB_DEPRECATED*` 与 `MATH_LIB_EXPORT`/`MATH_LIB_NO_EXPORT`/`MATH_LIB_DEPRECATED*` 是 DLL 构建宏，不形成 Rust ABI。
- `PYBIND11_MODULE(surfepy,m)` 是适配入口；`main` 是可视化 test harness 入口。

## 第三方与边界依赖图

| 依赖 | 冻结使用点 | 判定 |
| --- | --- | --- |
| C++11/STL | 全核心的 vector/set/sort/exception/iostream | 用 Rust 标准库与类型化错误替代 |
| Eigen Core/Dense/LU/Eigenvalues | 动态矩阵/向量、partial-pivot LU、LLT、特征分解、表达式运算 | 用纯 Rust 数组/稠密矩阵/求解和 3×3 特征路径替代；禁止 native BLAS/LAPACK |
| OpenMP | API 批量标量评估、最近邻统计、Greedy 残差与候选段 | 用安全 Rust 并发或确定性串行路径替代；不得保留数据竞争 |
| CMake | 构建 `math_lib`、`surfe_lib`、pybind module 和可选 geo_builder | 不迁移；T05 建 Cargo 护栏 |
| pybind11 | `surfepy` 绑定 | 不进入 Rust 生产核心；只保留 API 可达性证据 |
| Qt/VTK/geo_builder | `GEO_BUILDER` 可选分支和 `test/main.cpp` | 全部排除 |
| Windows DLL/console API | module headers、`debug.cpp` | ABI/调试排除 |

## 所有权摘要

- `Surfe_API` 独占一个裸指针 `method_`；具体模型独占或别名持有 `solver/kernel/rbf_kernel/p_basis` 裸指针。冻结源码未集中实现这些对象的完整析构所有权；Rust 目标必须用值、枚举或受控智能指针表达，而不复制泄漏/悬垂行为。
- `GRBF_Modelling_Methods::constraints` 是原始输入与后续原地更新的唯一模型容器；`interface_point_lists`、`interface_test_points`、`_increment_pairs` 是派生副本。
- 每个 solver 构造时复制矩阵/RHS，最终拥有 `weights`；评估读取模型约束、核和 `solver->weights`，把结果写回调用方拥有的 `Point`。
- `Modified_Kernel` 持有 RBF kernel 与 Lagrangian basis 指针，其复制/析构关系按源码记录；安全所有权重构属于实现任务，不在 T01 修改算法。

## 已登记的后续核验

- T02：对本清单每个符号做唯一能力分类，包括声明/定义差异、常量体、抛出体、TODO 与不可达性。
- T03/T04：以本清单和调用图定义 oracle/fixture 字段；本任务未构建 oracle 或生成 golden。
- T05：将 Eigen/OpenMP/CMake/pybind/native 边界落实为 Cargo 机器护栏。
- T06–T31：每个符号只映射到上述既有任务；未创建新任务。
