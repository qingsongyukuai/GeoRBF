# T02 冻结 Surfe 逐符号能力分类

## 身份、范围与分类规则

- 权威源码是 `https://github.com/MichaelHillier/surfe.git` 的提交
  `290dbe0ab344f4258a4935f05cad0f153f0f69a4`。
- 已逐一核对
  `surfe_lib/*.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、
  `math_lib/*.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、
  `surfe_pybindings/pybindings.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
  和 `test/main.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的声明、定义和活动调用点。
- “可达”只指从冻结的 `Surfe_API` 构造/公共方法或冻结 `test/main.cpp` 的活动核心调用
  出发可到达的非可视化路径。注释掉的调用、仅能直接实例化内部类后手工调用的成员、
  声明但无定义的成员，都不能把能力提升为可达。
- 分类码是互斥的。优先级依次为：范围决定 `V/X`；明确占位或哨兵体 `T`；公开可达
  的错误目标、未初始化读取、越界、漏返回或数据竞争 `D`；有有效子集但明确缺口 `P`；
  有定义却没有上述根到达的路径 `U`；其余为 `I`。

| 码 | 唯一分类 | 判定 |
| --- | --- | --- |
| `I` | 已实现 | 有具体定义，并在实际可达路径上产生所声明的非可视化行为。 |
| `P` | 部分实现 | 可达且有有效行为，但同一符号的声明/分支包含明确缺口。 |
| `T` | TODO | 空体、恒定成功占位、`TO IMPLEMENT`，或 `throw -666` 哨兵体。 |
| `U` | 不可达 | 有声明或定义，但冻结活动调用图没有从上述根到它的边。 |
| `D` | 缺陷 | 存在可定位的错误行为或 C++ 未定义行为；是否公开可达另由调用证据说明，只记录、不在 T02 修复。 |
| `V` | 可视化 | 控制台调试、VTK/GUI 显示或纯可视化 harness。 |
| `X` | 排除 | C++ ABI、构建/绑定适配或 `geo_builder` 边界，不形成 Rust 核心能力。 |

下文的 `{a,b}` 展开为两个独立符号；`C × M` 展开为每个 `C::m`。每个展开后的符号
恰好出现在一个分类行，并只给出一个“关闭任务”。字段也按符号分类；局部变量不属于
T01/T02 的公开声明集。

## 参数、枚举、常量和参数字段

枚举值集合精确为：

- `SD = {DXDX,DXDY,DXDZ,DYDX,DYDY,DYDZ,DZDX,DZDY,DZDZ}`；
- `FD = {DX,DY,DZ}`；
- `RK = {Cubic,Gaussian,MQ,MQ3,IMQ,TPS,R,WendlandC2,MaternC4}`；
- `ST = {Linear,Quadratic}`；
- `MT = {Single_surface,Lajaunie_approach,Stratigraphic_horizons,Continuous_property,Vector_field}`；
- `AV = {Xaxis,Yaxis,Zaxis}`。

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `D2R`、`R2D`、`Epilson` | `I` | 方向换算、角残差和 `collocated` 直接读取；见 `modelling_parameters.h`、`modelling_input.{h,cpp}@commit`。 | T06 |
| `Parameter_Types::{SecondDerivatives,FirstDerivatives,RBF,SolverType,ModelType}` 及 `SD`、`FD`、`RK`、`ST`、`MT` 的每个展开值 | `I` | 矩阵分量、核工厂、solver 和模型工厂均有活动分支。 | T06 |
| `Parameter_Types::DWRT`、`Parameter_Types::{PT1,PT2}` | `U` | 仅声明；全范围无读取点。 | T06 |
| `Parameter_Types::AXIS` 及 `AV` 的每个展开值 | `U` | 只被同样不可达的 `get_maximal_axial_variability_order` 使用。 | T09 |
| `Parameters::Parameters` | `I` | 五模型默认构造和 `Surfe_API(Parameters)` 使用。 | T06 |
| `Parameters::{model_type,min_stratigraphic_thickness,basis_type,shape_parameter,polynomial_order,model_global_anisotropy,use_restricted_range,smoothing_amount,use_regression_smoothing,interface_uncertainty,angular_uncertainty}` | `I` | 至少一个活动工厂、装配、restricted 或 smoothing 分支读取。 | T06 |
| `Parameters::{use_interface,use_planar,use_tangent,use_inequality}` | `U` | Add 方法只写入；冻结核心无读取点。 | T06 |
| `Parameters::advanced_parameters` | `U` | 只声明并默认初始化，无写入或读取点。 | T06 |
| `Parameters::use_greedy` | `U` | 默认值和 setter 写入后无任何消费点。 | T31 |
| `InternalParameters::InternalParameters` 及 `InternalParameters::{n_interface,n_planar,n_inequality,n_tangent,n_constraints,n_equality,modified_basis,poly_term,n_poly_terms,problem_type,restricted_range}` | `I` | 五模型 layout、basis、solver、装配和评估直接读写。 | T16 |
| `InputParameters` 及 `InputParameters::{parameters,interface_file,planar_file,tangent_file,inequality_file}` | `U` | 完整范围只有声明，无构造或读取点。 | T06 |

这里的 `@commit` 均指本文件首节固定的完整提交。`T01` 文本所说“22 个异常”是计数
笔误：该文件实际声明并列出了 23 个异常类型；没有漏掉异常符号。

## 约束值类型

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `Point::{_x,_y,_z,_c,_scalar_field,_field_normal}`、`Point::Point(x,y,z,c)`、`Point::{x,y,z,c,scalar_field,set_scalar_field,set_vector_field,nx_interp,ny_interp,nz_interp,operator<}`、`collocated` | `I` | API 输入、排序、核半径和评估返回直接使用。 | T07 |
| `Point::Point()` | `U` | 空默认构造体存在但活动核心不调用；其字段未初始化，Rust 不复制该无效输入行为。 | T07 |
| `Point::{set_x,set_y,set_z}` | `I` | Lagrangian 轴平面回退可达调用。 | T11 |
| `Point::set_c` | `U` | 只有声明/内联定义，无调用点。 | T07 |
| `Interface::{_level,_level_bound,Interface(x,y,z,level,c),level,level_lower_bound,level_upper_bound,setLevelBounds}` | `I` | 分组、RHS 和 restricted bounds 可达。 | T07 |
| `Interface::{_residual,residual,setResidual}` | `U` | 只在不可达 Greedy 链中读取/写入。 | T31 |
| `Interface::setLevel` | `U` | 只在无调用点的 reconstruction 函数体中使用。 | T21 |
| `Interface::Interface()` | `U` | 空默认构造无活动调用点。 | T07 |
| `Inequality::{_inequality_level,Inequality(x,y,z,level,c),level}` | `I` | Single/Stratigraphic 普通和 restricted 分支读取。 | T07 |
| `Inequality::{_residual,residual,setResidual}` | `U` | 只在不可达 Greedy 链中使用。 | T31 |
| `Inequality::Inequality()` | `U` | 空默认构造无活动调用点。 | T07 |
| `Planar::{_dip,_strike,_polarity,_normal,_normal_bound,Planar(normal),Planar(dip,strike,polarity),_compute_strike_dip_polarity_from_normal,_compute_normal_from_strike_dip_polarity,nx,ny,nz,nx_lower_bound,nx_upper_bound,ny_lower_bound,ny_upper_bound,nz_lower_bound,nz_upper_bound,setNormalBounds}` | `I` | 两种 API 构造、各模型 RHS 和 restricted bounds 可达。 | T07 |
| `Planar::{dip,strike,polarity,getDipVector,getStrikeVector}` | `U` | 有定义但活动核心不调用这些访问器。 | T07 |
| `Planar::{_residual,residual,setResidual}` | `U` | 只服务不可达 Greedy。 | T31 |
| `Planar::setNormal` | `U` | 只服务不可达 reconstruction。 | T21 |
| `Planar::Planar()` | `U` | 空默认构造无活动调用点。 | T07 |
| `Tangent::{_tangent,_angle_bound,_inner_product_constraint,Tangent(x,y,z,t),tx,ty,tz,angle_lower_bound,angle_upper_bound,inner_product_constraint,setAngleBounds}` | `I` | Tangent 装配、评估和 restricted bounds 可达。 | T07 |
| `Tangent::{_residual,residual,setResidual}` | `U` | 只服务不可达 Greedy。 | T31 |
| `Tangent::setInnerProductConstraint` | `U` | 只服务不可达 reconstruction。 | T21 |
| `Tangent::Tangent()` | `U` | 空默认构造无活动调用点。 | T07 |
| `Constraints::{inequality,itrface,planar,tangent,Constraints}` | `I` | API 与全部模型的输入容器。 | T07 |
| `Constraints::{_avg_nn_dist_ie,_avg_nn_dist_itr,_avg_nn_dist_p,_avg_nn_dist_t,compute_inequality_avg_nn_distance,compute_interface_avg_nn_distance,compute_planar_avg_nn_distance,compute_tangent_avg_nn_distance,compute_avg_nn_distances,GetInequalityAvgNNDist,GetInterfaceAvgNNDist,GetPlanarAvgNNDist,GetTangentAvgNNDist,SetInequalityAvgNNDist,SetInterfaceAvgNNDist,SetPlanarAvgNNDist,SetTangentAvgNNDist}` | `U` | 唯一调用者是不可达 `run_greedy_algorithm` 及其 hooks。 | T31 |
| `SpatialParameters::{resolution,xmin,xmax,ymin,ymax,zmin,zmax}` | `I` | `GetDataBoundsAndResolution` 填充并返回。 | T09 |

`Point::operator==` 在 `PLAN.md` 的 T08 来源提示中出现，但冻结头、定义和调用点中均不存在；
它不是待迁移符号。实际同点语义只有自由函数 `collocated`。

## 空间和 Greedy 候选辅助

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `convert_constraints_to_points`、`distance_btw_pts`、`get_largest_distance_between_points` | `I` | Single/Stratigraphic restricted range 和公共空间查询可达。 | T09 |
| `spatial_metrics` | `D` | 空输入仍返回 `true`，把 bounds 留为 `±DBL_MAX`；所以公共 `GetDataBoundsAndResolution` 不会抛其声明的空间错误。 | T09 |
| `nearest_neighbour_index`、`get_n_nearest_neighbours_to_point`、两个 `furtherest_neighbour_index`、`avg_nn_distance`、`Find_STL_Vector_Indices_FurtherestTwoPoints`、`Find_STL_Vector_Index_ofPointClosestToOtherPointWithinDistance`、`calculate_bounds`、`get_extremal_point_data_indices_from_points`、`is_index_in_list`、`get_maximal_axial_variability_order` | `U` | 有具体定义，但只被不可达 Greedy 链调用或完全无调用点。 | T09 |
| `Get_Inequality_STL_Vector_Indices_With_Large_Residuals`、`Get_Interface_STL_Vector_Indices_With_Large_Residuals`、`Get_Planar_STL_Vector_Indices_With_Large_Residuals` | `U` | 只被不可达模型 Greedy hooks 调用。 | T31 |
| `Get_Tangent_STL_Vector_Indices_With_Large_Residuals` | `D` | 除 Greedy 不可达外，存在“大残差集合为空”时走到函数末尾而不返回值的确定漏返回。 | T31 |

## Kernel、RBF、Polynomial 与 Modified Kernel

定义方法集：

- `K = {basis,dx_p1,dx_p2,dy_p1,dy_p2,dz_p1,dz_p2,dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz,clone}`。
- `KD = K - {basis,clone}`。
- `CI = {Cubic,Gaussian,MQ,MQ3,TPS,IMQ,WendlandC2,MaternC4}`。
- `CA = {ACubic,AGaussian,AMQ,ATPS,AIMQ}`。

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `Kernel::{_p1,_p2,set_points,p1,p2,basis_pt_pt,basis_pt_planar_x,basis_planar_x_pt,basis_pt_planar_y,basis_planar_y_pt,basis_pt_planar_z,basis_planar_z_pt,basis_pt_tangent,basis_tangent_pt,basis_planar_planar,basis_tangent_tangent,basis_planar_tangent,basis_tangent_planar,clone,~Kernel}` | `I` | 抽象派发槽均有 RBF/Modified 覆写；装配和评估可达。 | T15 |
| `Kernel::{set_center,set_evaluation_points}` | `U` | 两个内联 setter 有定义，但全范围无调用点。 | T15 |
| `RBFKernel::{_radius,_x_delta,_y_delta,_z_delta,_c_delta,_Transform,RBFKernel,~RBFKernel,radius,scaled_radius,get_global_anisotropy,basis,basis_pt_pt,basis_pt_planar_x,basis_planar_x_pt,basis_pt_planar_y,basis_planar_y_pt,basis_pt_planar_z,basis_planar_z_pt,basis_pt_tangent,basis_tangent_pt,basis_planar_planar,basis_tangent_tangent,basis_planar_tangent,basis_tangent_planar,clone}` | `I` | 普通/anisotropic 工厂与模型矩阵可达。 | T12 |
| `RBFKernel::{dx_p1,dx_p2,dy_p1,dy_p2,dz_p1,dz_p2,dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz}` | `P` | 抽象槽对大多数核完整，但 `R`/`AR` 的对应覆写全是哨兵体。 | T12 |
| `RBFKernel::_Global_Plunge` | `U` | anisotropy 写入后没有任何读取点。 | T13 |
| `CI × K` | `I` | 每个展开符号均有具体定义；`MQ3`、`WendlandC2`、`MaternC4` 也各有 17 个定义，并被普通核工厂选择。 | T12 |
| `CA × K` | `I` | 每个展开符号均有具体 anisotropic 定义并被工厂选择。 | T13 |
| `R::{basis,clone,~R}` | `I` | `basis` 返回半径，普通核工厂可选择。 | T12 |
| `R × KD` | `T` | 15 个展开符号中 12 个直接 `throw -666`，3 个别名 Hessian 转入这些哨兵。 | T12 |
| `AR::{AR,basis,clone,~AR}` | `I` | anisotropic 工厂可选择并计算 scaled radius。 | T13 |
| `AR × KD` | `T` | 15 个展开符号直接或间接进入 `throw -666`。 | T13 |
| `Gaussian::_shape_parameter`、`MQ::_shape_parameter`、`MQ3::_shape_parameter`、`IMQ::_shape_parameter`、`AGaussian::_shape_parameter`、`AMQ::_shape_parameter`、`AIMQ::_shape_parameter`、`WendlandC2::_cutoff`、`MaternC4::_s` | `I` | 冻结声明集中的私有核参数，构造器和核体读取。 | T12 |
| `Cubic::~Cubic`、`Gaussian::{Gaussian,~Gaussian}`、`MQ::{MQ,~MQ}`、`MQ3::{MQ3,~MQ3}`、`TPS::~TPS`、`IMQ::{IMQ,~IMQ}`、`WendlandC2::{WendlandC2,~WendlandC2}`、`MaternC4::{MaternC4,~MaternC4}` | `I` | 与 `CI × K` 对应的显式生命周期符号。 | T12 |
| `ACubic::{ACubic,~ACubic}`、`AGaussian::{AGaussian,~AGaussian}`、`AMQ::{AMQ,~AMQ}`、`ATPS::{ATPS,~ATPS}`、`AIMQ::{AIMQ,~AIMQ}` | `I` | 与 `CA × K` 对应的显式生命周期符号。 | T13 |
| `Polynomial_Basis::{_p,_truncated,Polynomial_Basis,set_point,basis,dx,dy,dz,clone}` | `I` | 三个具体 polynomial 完整覆写，普通模型装配/评估可达。 | T10 |
| `Poly_Zero::{Poly_Zero,basis,dx,dy,dz,clone}`、`Poly_First::{Poly_First,basis,dx,dy,dz,clone}`、`Poly_Second::{Poly_Second,basis,dx,dy,dz,clone}` | `I` | 普通与 truncated 0/1/2 阶都有具体项序定义。 | T10 |
| `Lagrangian_Polynomial_Basis::{_polynomial_constants,_derivative_polynomial_constants,unisolvent_subset_points,poly,poly_dx,poly_dy,poly_dz}` | `I` | Modified Kernel 构造和所有组合读取。 | T11 |
| `Lagrangian_Polynomial_Basis::{Lagrangian_Polynomial_Basis,_get_unisolvent_subset,_initialize_basis}` | `P` | 只处理一阶且只识别轴平面退化；相等最大轴会执行多个独立 `if` 并可能选出多于四点，任意共面 determinant 也未验证。 | T11 |
| `Modified_Kernel::{_aRBFKernel,_aLPB,Modified_Kernel,Modified_Kernel(copy),~Modified_Kernel,basis_pt_pt,basis_pt_planar_x,basis_planar_x_pt,basis_pt_planar_y,basis_planar_y_pt,basis_pt_planar_z,basis_planar_z_pt,basis_pt_tangent,basis_tangent_pt,basis_planar_planar,basis_tangent_tangent,basis_planar_tangent,basis_tangent_planar,clone}` | `I` | 普通 QP/restricted 装配和评估可达；全部 Value/Planar/Tangent 组合有定义。 | T14 |

`CI × K` 展开为 136 个符号，`CA × K` 展开为 85 个符号，`R × KD` 和
`AR × KD` 各展开为 15 个符号。这四行没有交集。源码扫描找到 `R/AR` 的 24 个
字面 `throw -666` 体和 6 个转入它们的别名 Hessian，合计 30 个哨兵导数符号；
`MQ3`、`WendlandC2`、`MaternC4` 不是 TODO 或恒定返回实现。

## 共享模型基类

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `GRBF_Modelling_Methods::{intern_params,interface_iso_values,interface_point_lists,constraints,parameters,solver,kernel,rbf_kernel,interface_test_points}` | `I` | 五模型的输入、派生 layout、basis、solver 和输出状态。 | T16 |
| `GRBF_Modelling_Methods::{_iteration,error_msg}` | `U` | `_iteration` 只服务不可达 Greedy；`error_msg` 无读取/写入点。 | T31 |
| `GRBF_Modelling_Methods::{_get_distinct_interface_iso_values,_get_interface_points,_get_distinct_inequality_iso_values,get_interface_data,_update_interface_iso_values,remove_collocated_constraints,setup_basis_functions,get_equality_matrix,create_rbf_kernel}` | `I` | `ComputeInterpolant` 或具体模型活动路径直接调用。 | T16 |
| `GRBF_Modelling_Methods::check_input_data` | `P` | 活动 Stratigraphic 路径可达，但四类检查中只有 inequality/interface level 冲突有实现。 | T30 |
| `GRBF_Modelling_Methods::_interface_points_are_coplanar` | `T` | 内联恒定 `true`，并明确注释 “Not implemented yet”。 | T11 |
| `GRBF_Modelling_Methods::{get_interface_points_ouput,get_interface_iso_values,check_interpolant,get_method,run_greedy_algorithm,_SetIteration,~GRBF_Modelling_Methods}` | `U` | 无冻结根调用；`get_method` 只被同样不可达的 `run_greedy_algorithm` 调用。 | T31 |
| `GRBF_Modelling_Methods::_output_greedy_debug_objects` | `V` | 唯一调用被注释，函数体只有调试输出占位检查。 | T31 |
| 基类虚槽 `GRBF_Modelling_Methods::{get_interpolation_matrix,get_equality_values,eval_scalar_interpolant_at_point,eval_vector_interpolant_at_point,get_method_parameters,setup_system_solver}` | `I` | `ComputeInterpolant`/Evaluate 经五模型虚派发到具体定义。 | T16 |
| 基类虚槽 `GRBF_Modelling_Methods::process_input_data` | `P` | 四模型有工作体，Vector Field 覆写为空。 | T16 |
| 基类虚槽 `GRBF_Modelling_Methods::{get_minimial_and_excluded_input,measure_residuals,append_greedy_input,convert_modified_kernel_to_rbf_kernel,clone}` | `U` | Greedy、reconstruction 和模型 clone 均无冻结根调用。 | T31 |

## 五个具体模型

### Single Surface

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `Single_Surface::{Single_Surface(),Single_Surface(Parameters),_get_polynomial_matrix_block,_insert_polynomial_matrix_blocks_in_interpolation_matrix,create_polynomial_basis,get_equality_values,process_input_data}`、`Single_Surface::p_basis` | `I` | 普通 equality/poly 路径可达。 | T22 |
| `Single_Surface::{get_interpolation_matrix,setup_system_solver,eval_scalar_interpolant_at_point}` | `I` | Linear、ordinary QP 和 LOQO 分支均有活动调用。 | T24 |
| `Single_Surface::get_method_parameters` | `D` | `intern_params.restricted_range` 只置 `true` 不复位；同一对象从 restricted 切回 false 后仍走 LOQO。 | T24 |
| `Single_Surface::eval_vector_interpolant_at_point` | `D` | 创建 `kernel_j` 却全部用共享 `kernel` 求值；并发评估会竞争 `_p1/_p2`。 | T29 |
| `Single_Surface::{get_inequality_matrix,get_inequality_values(vector)}` | `P` | ordinary inequality 子路径可达；各自的 restricted 子分支没有调用点。 | T23 |
| `Single_Surface::get_inequality_values(b,r)` | `I` | restricted LOQO 活动路径调用。 | T24 |
| `Single_Surface::{get_minimial_and_excluded_input,measure_residuals,append_greedy_input}` | `U` | 只有不可达 `run_greedy_algorithm` 能调用。 | T31 |
| `Single_Surface::convert_modified_kernel_to_rbf_kernel` | `U` | 有完整定义但全仓库零调用点。 | T21 |
| `Single_Surface::{clone,~Single_Surface}` | `U` | 模型 clone/析构没有活动调用；`Surfe_API` 也未声明析构。 | T29 |

### Lajaunie Approach

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `Lajaunie_Approach::{Lajaunie_Approach(),Lajaunie_Approach(Parameters),_get_increment_pairs,_get_polynomial_matrix_block,_insert_polynomial_matrix_blocks_in_interpolation_matrix,create_polynomial_basis,get_equality_values,get_inequality_values(b,r),get_interpolation_matrix,process_input_data,eval_scalar_interpolant_at_point}`、`Lajaunie_Approach::{_n_increment_pair,_increment_pairs,p_basis}` | `I` | 普通和 restricted Lajaunie 活动路径调用。 | T25 |
| `Lajaunie_Approach::get_method_parameters` | `D` | restricted 状态只置位不复位，与 Single 相同。 | T25 |
| `Lajaunie_Approach::setup_system_solver` | `D` | LOQO 失败错误地抛 `pc_quadratic_solver_failure`，不是已声明的 LOQO 类别。 | T25 |
| `Lajaunie_Approach::eval_vector_interpolant_at_point` | `D` | clone 后仍使用共享 `kernel`，存在 `_p1/_p2` 数据竞争。 | T29 |
| `Lajaunie_Approach::{get_minimial_and_excluded_input,measure_residuals,append_greedy_input}` | `U` | 只有不可达 Greedy 调用。 | T31 |
| `Lajaunie_Approach::convert_modified_kernel_to_rbf_kernel` | `U` | 有完整定义但零调用点。 | T21 |
| `Lajaunie_Approach::{clone,~Lajaunie_Approach}` | `U` | 无活动调用。 | T29 |

### Stratigraphic Surfaces

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `Stratigraphic_Surfaces::{Stratigraphic_Surfaces(),Stratigraphic_Surfaces(Parameters),_get_lithostratigraphic_increment_pairs_for_inequality_point,_get_closest_horizon_level_above_given_level,_get_closest_horizon_level_below_given_level,get_equality_values,get_inequality_values(vector),get_inequality_values(b,r),get_inequality_matrix,get_interpolation_matrix,setup_system_solver,eval_scalar_interpolant_at_point}`、`Stratigraphic_Surfaces::{_n_increment_pairs,_n_sequenced_interface_pairs,_n_sequenced_inequality_pairs,_n_interface_pairs}` | `I` | ordinary QP 和 restricted LOQO 活动路径调用。 | T26 |
| `Stratigraphic_Surfaces::_get_increment_pairs`、`Stratigraphic_Surfaces::process_input_data`、`Stratigraphic_Surfaces::_increment_pairs` | `D` | `_increment_pairs` 在重复拟合前不清空；`resize` 保留旧 pair 后再次 `push_back`，且零 pair 仍返回成功。 | T26 |
| `Stratigraphic_Surfaces::get_method_parameters` | `D` | restricted 状态只置位不复位。 | T26 |
| `Stratigraphic_Surfaces::eval_vector_interpolant_at_point` | `D` | clone 后仍使用共享 `kernel`，存在 `_p1/_p2` 数据竞争。 | T29 |
| `Stratigraphic_Surfaces::{_get_polynomial_matrix_block,_insert_polynomial_matrix_blocks_in_interpolation_matrix,create_polynomial_basis,p_basis}` | `U` | `poly_term` 只在不可达 reconstruction 体中开启。 | T21 |
| `Stratigraphic_Surfaces::{get_minimial_and_excluded_input,measure_residuals,append_greedy_input}` | `T` | 三个内联体均恒定返回 `true` 且标注 `TO implement`。 | T31 |
| `Stratigraphic_Surfaces::convert_modified_kernel_to_rbf_kernel` | `U` | 有完整定义但零调用点。 | T21 |
| `Stratigraphic_Surfaces::{clone,~Stratigraphic_Surfaces}` | `U` | 无活动调用。 | T29 |

### Continuous Property

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `Continuous_Property::{Continuous_Property(),Continuous_Property(Parameters),process_input_data,setup_system_solver,eval_scalar_interpolant_at_point,get_interpolation_matrix}` | `I` | 公开工厂可达的 interface-value LU 路径。 | T27 |
| `Continuous_Property::get_method_parameters` | `P` | 明确把 planar/tangent/inequality/poly 计数置零，只实现 interface-value 子集。 | T27 |
| `Continuous_Property::get_equality_values` | `D` | 输出向量按 interface 计数分配，却遍历实际 planar/tangent 容器；公共 API 添加这些约束时可越界写。 | T27 |
| `Continuous_Property::eval_vector_interpolant_at_point` | `D` | clone 后仍使用共享 `kernel`，存在 `_p1/_p2` 数据竞争。 | T29 |
| `Continuous_Property::{_get_polynomial_matrix_block,_insert_polynomial_matrix_blocks_in_interpolation_matrix,create_polynomial_basis,p_basis}` | `U` | 活动参数路径恒令 `n_poly_terms=0, poly_term=false`。 | T27 |
| `Continuous_Property::get_minimial_and_excluded_input` | `T` | 恒定返回 `true`，相邻注释保留了未完成的 intended body。 | T31 |
| `Continuous_Property::{measure_residuals,append_greedy_input}` | `U` | 只有不可达 Greedy 调用。 | T31 |
| `Continuous_Property::convert_modified_kernel_to_rbf_kernel` | `T` | 恒定返回 `true` 并标注 `To IMPLEMENT`。 | T27 |
| `Continuous_Property::{clone,~Continuous_Property}` | `U` | 无活动 clone/delete 调用；析构中的 `dest` 输出不构成可达能力。 | T29 |

### Vector Field

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `Vector_Field::{Vector_Field(),Vector_Field(Parameters),get_method_parameters,get_equality_values,get_interpolation_matrix,setup_system_solver,eval_scalar_interpolant_at_point,eval_vector_interpolant_at_point}` | `I` | 公开五分支工厂可达的 Planar Hessian、势和梯度路径。 | T28 |
| `Vector_Field::process_input_data` | `P` | 公开可达但为空体；没有 `no_planar_data` 检查，零 Planar 会继续进入零维求解。 | T28 |
| `Vector_Field::{get_minimial_and_excluded_input,measure_residuals,append_greedy_input,convert_modified_kernel_to_rbf_kernel}` | `T` | 四个内联体均恒定 `true` 并标注 `TO implement`。 | T31 |
| `Vector_Field::{clone,~Vector_Field}` | `U` | 无活动调用。 | T29 |

五模型的 12 个覆写槽均已在以上各表出现一次。尤其是完整 reconstruction 函数体不能因
“有定义”升级为能力：`single_surface.cpp`、`lajaunie.cpp`、`stratigraphic_surfaces.cpp@commit`
中的三个 `convert_modified_kernel_to_rbf_kernel` 在声明/定义/调用集合核对中调用集合为空。

## Solver 与数学辅助

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `System_Solver::{weights,System_Solver,solve}` | `I` | 三个 concrete solver 覆写，模型拟合读取 weights。 | T18 |
| `System_Solver::{~System_Solver,validate_matrix_systems}` | `U` | solver 裸指针不 delete；三个 validate 覆写也没有活动调用。 | T18 |
| `Linear_LU_decomposition::{_interpolation_matrix,_constraint_values,Linear_LU_decomposition(matrix,vector),solve}` | `I` | 四模型 Linear 路径调用。 | T18 |
| `Linear_LU_decomposition::{Linear_LU_decomposition(),~Linear_LU_decomposition,validate_matrix_systems,check_solution}` | `U` | 默认构造、析构和两个诊断函数均无活动调用。 | T18 |
| `Quadratic_Predictor_Corrector::{_interpolation_matrixD,_hessian_matrixD,_equality_matrixD,_inequality_matrixD,_equality_vectorD,_inequality_vectorD,Quadratic_Predictor_Corrector,solve}` | `I` | Single/Stratigraphic ordinary QP 调用。 | T19 |
| `Quadratic_Predictor_Corrector::validate_matrix_systems` | `U` | `solve` 中调用被注释。 | T19 |
| `Quadratic_Predictor_Corrector_LOQO::{_H,_A,_b,_r,Quadratic_Predictor_Corrector_LOQO,solve}` | `I` | 三模型 restricted 路径调用。 | T20 |
| `Quadratic_Predictor_Corrector_LOQO::validate_matrix_systems` | `U` | 无调用点。 | T20 |
| `Math_methods::{_find_step_length,max_element_wrt_zero,quadratic_solver}` | `I` | ordinary predictor-corrector 直接调用。 | T19 |
| `Math_methods::{_find_step,_find_positivity_step,quadratic_solver_loqo}` | `I` | LOQO 路径直接调用。 | T20 |
| `Math_methods::{sort_vector_w_index,SWAP}` | `I` | 分组、Lagrangian 和地层邻层选择可达。 | T08 |
| `Math_methods::_get_double` | `I` | ordinary QP 活动体调用；只影响冻结进度诊断值。 | T19 |
| `Math_methods::{_rot,RandomDouble}` | `U` | `_rot` 只有声明且无定义；`RandomDouble` 有定义但零调用。 | T19 |
| `Math_methods::angle_btw_2_vectors` | `U` | 唯一调用者是不可达 Greedy residual hooks。 | T31 |

## `Surfe_API` 状态机和公共入口

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `Surfe_API::method_`、`Surfe_API::get_method_from_parameters`、`Surfe_API::Surfe_API(Parameters)` | `I` | 五模型公开工厂和安全初始化路径完整。 | T29 |
| `Surfe_API::{have_interpolant_,parameters_changed_,constraints_changed_,Surfe_API(int),InterpolantComputed}` | `D` | int 构造只写 `parameters_changed_`，其余两个 bool 未初始化；`InterpolantComputed` 可直接读取不确定值。 | T30 |
| `Surfe_API::{AddInterfaceConstraint,AddPlanarConstraintwNormal,AddPlanarConstraintwStrikeDipPolarity,AddPlanarConstraintwAzimuthDipPolarity,AddTangentConstraint,AddInequalityConstraint}` | `I` | 参数化约束构造并写相应容器。 | T29 |
| `Surfe_API::ComputeInterpolant` | `P` | 正常 process/basis/solver 路径完整，但从不读取 `use_greedy`，也不调用 reconstruction。 | T29 |
| `Surfe_API::SetRegressionSmoothing` | `D` | 忽略传入 bool，无条件把 `use_regression_smoothing` 置 `true`。 | T30 |
| `Surfe_API::SetGreedyAlgorithm` | `D` | 忽略传入 bool、无条件置 `true`，且该字段没有消费点。 | T31 |
| `Surfe_API::{SetRestrictedRange,SetRBFKernel(enum),SetRBFKernel(name),SetRBFShapeParameter,SetPolynomialOrder,SetGlobalAnisotropy}` | `I` | 均有参数写入和后续活动读取或精确拒绝分支。 | T30 |
| `Surfe_API::{EvaluateInterpolantAtPoint,EvaluateVectorInterpolantAtPoint}` | `P` | Parameters 构造路径有完整状态检查；int 构造的未初始化 bool 使拟合前行为不可靠。 | T30 |
| `Surfe_API::{EvaluateInterpolantAtPoints,EvaluateVectorInterpolantAtPoints}` | `D` | 与单点不同，拟合后参数/约束变化时不检查 changed flags；scalar 批量进度计数另有 OpenMP 数据竞争。 | T30 |
| `Surfe_API::GetDataBoundsAndResolution` | `D` | 依赖缺陷 `spatial_metrics`，空约束不会得到预期错误。 | T30 |
| `Surfe_API::{GetInterfaceReferencePoints,GetInterfaceConstraints,SetInterfaceConstraints,GetPlanarConstraints,SetPlanarConstraints,GetTangentConstraints,GetInequalityConstraints,SetInequalityConstraints,GetNumberOfInterfaces}` | `I` | 定义完整并从公共 API/pybind 可达。 | T29 |
| `Surfe_API::SetTangentConstraints` | `D` | 清空 tangent 后逐行调用 `AddPlanarConstraintwNormal`，数据进入 planar 而非 tangent。 | T30 |
| `Surfe_API::progress` | `V` | 只输出控制台百分比；不属于数值/API 结果。 | T29 |

`Surfe_API` 没有析构声明，因此 `method_` 及模型持有的 solver/kernel/polynomial 裸指针在正常
API 生命周期中不释放。这是 C++ 所有权缺陷，不是一个可分类的已声明符号；T29 必须以 Rust
所有权安全替代，不能复制泄漏。

## 异常、调试、构建和适配边界

定义已使用异常集合：

`EU = {nointerfacedata,nointerfaceincrementpairs,invalidinputdata,
failurecomputingglobalanisotropy,failurecreatinganisotropickernel,
failuresettingupbasisfunctions,failurecreatingmodifiedkernel,
failurecreatinglagrangianpolynomialbasis,linearsolverfailure,
pcquadratricsolverfailure,loqoquadratricsolverfailure,
errorcomputinginterpolationmatrix,errorcomputingequalityvector,
errorcomputinginequalityvector,errorupdatinginterfaceisovalues,
missinginterpolant,unknownrbf,interpolantneedsupdate,unknownmodellingmode,
problemcomputingspatialparameters,arrayhasincorrectdimensions}`。

它们对应的已使用常量集合为：

`EC = {no_iterface_data,no_interface_increment_pairs,invalid_input_data,
failure_computing_global_anisotropy,failure_creating_anisotropic_kernel,
failure_setting_up_basis_functions,failure_creating_modified_kernel,
failure_creating_lagrangian_polynomial_basis,linear_solver_failure,
pc_quadratic_solver_failure,loqo_quadratic_solver_failure,
error_computing_interpolation_matrix,error_computing_equality_vector,
error_computing_inequality_vector,error_updating_interface_iso_values,
missing_interpolant,unknown_rbf,interpolant_needs_update,unknown_modelling_mode,
problem_computing_spatial_parameters,array_has_incorrect_dimensions}`。

| 符号集合 | 分类 | 定义/调用证据 | 关闭任务 |
| --- | --- | --- | --- |
| `EU` 中每个异常类型及其 `what`，以及 `GRBF_Exceptions::EC` 的每个展开常量 | `I` | 至少一个活动 throw/catch 路径引用。 | T06 |
| `noplanardata`、`noplanardata::what`、`GRBF_Exceptions::no_planar_data` | `U` | Vector Field 空 `process_input_data` 未使用该已声明错误。 | T06 |
| `errorcomputinginterpolant`、`errorcomputinginterpolant::what`、`GRBF_Exceptions::error_computing_interpolant` | `U` | 全范围无 throw 点。 | T06 |
| `SurfeExceptions::{errors,append_exceptions,SurfeExceptions,what}` | `I` | `ComputeInterpolant` 三段 catch 负责扁平化 nested exceptions。 | T06 |
| `open_console_window` | `V` | `debug.cpp@commit` 仅 Windows console 分配；无核心调用点。 | 明确排除 |
| `SURFE_LIB_EXPORT`、`SURFE_LIB_NO_EXPORT`、`SURFE_LIB_DEPRECATED`、`SURFE_LIB_DEPRECATED_EXPORT`、`SURFE_LIB_DEPRECATED_NO_EXPORT`、`SURFE_LIB_NO_DEPRECATED`、`DEFINE_NO_DEPRECATED` | `X` | `surfe_lib_module.h@commit` 的 DLL ABI 宏。 | 明确排除 |
| `MATH_LIB_EXPORT`、`MATH_LIB_NO_EXPORT`、`MATH_LIB_DEPRECATED`、`MATH_LIB_DEPRECATED_EXPORT`、`MATH_LIB_DEPRECATED_NO_EXPORT`、`MATH_LIB_NO_DEPRECATED` | `X` | `math_lib_module.h@commit` 的 DLL ABI 宏。 | 明确排除 |
| `GEO_BUILDER`、`_WIN32` 条件宏 | `X` | 分别只选择排除的 geo_builder 构建路径和 Windows console 适配。 | 明确排除 |
| `PYBIND11_MODULE(surfepy,m)` | `X` | 只把 30 个已有 `Surfe_API` 入口适配到 Python；pybind11 不迁移。 | 明确排除 |
| `main` | `V` | 活动核心前缀后只构建 grid、写 VTK 并显示；harness 不迁移。 | 明确排除 |
| `geo_builder/**` 的 `Geo_Builder` 等全部类型/函数、Qt/VTK/GUI、regular-grid/isosurface 写出与显示 | `X` | 不在 T01 非可视化声明集，仅作为 `test/main.cpp` 的排除边界。 | 明确排除 |

`debug.h` 中整段 `outc` templates 和 `matrix_solver.cpp` 中 GMP 转换体均被注释，不形成 C++
声明/定义符号；不能分类成能力。

## 缺陷和不可达性的可观察含义

T02 尚未建立 T03 oracle，因此下表是可机器复核的冻结源码证据，不声称数值 oracle 已验证。
后续任务必须先保存 oracle/兼容用例，再决定 Rust 的安全错误；不得复制 UB 或数据竞争。

| 证据 | 冻结可观察含义 | 唯一后续归属 |
| --- | --- | --- |
| `surfe_api.cpp::SetTangentConstraints@commit` | 合法 `n×6` tangent 矩阵清空旧 tangent，却增加 planar；随后 `GetTangentConstraints` 为空而 `GetPlanarConstraints` 增长。 | T30 |
| `surfe_api.cpp::Surfe_API(int)@commit` | `InterpolantComputed` 或拟合前 Evaluate 读取未初始化状态，结果不稳定；Rust 必须稳定拒绝拟合前评估。 | T30 |
| `surfe_api.cpp::{SetRegressionSmoothing,SetGreedyAlgorithm}@commit` | 传 `false` 仍启用相应 flag；其中 Greedy flag 再无读取点。 | T30 |
| `surfe_api.cpp::Evaluate*AtPoints@commit` | 拟合后修改参数/约束仍走旧 weights；标量 OpenMP 进度状态还有数据竞争。 | T30 |
| `modeling_methods.cpp::run_greedy_algorithm@commit` 和全仓调用集合 | 有完整循环体但没有 API 调用边，且即使直接调用也不把 `greedy_method` 结果装回当前对象。不能宣称 Greedy 可用。 | T31 |
| 三个 `convert_modified_kernel_to_rbf_kernel@commit` 和全仓调用集合 | Single/Lajaunie/Stratigraphic 定义完整但调用集合为空；冻结公开结果是 QP/LOQO weights 直接评估，不是二次 LU reconstruction。 | T21 |
| `basis.cpp::R/AR::{dx_*,d*,...}@commit` | `R`/`AR` 可做 point/point value；任何需要一阶或混合二阶导数的可达约束/梯度路径抛出非 `std::exception` 的整数 `-666`。 | T12 |
| `modelling_input.cpp::Get_Tangent_STL_Vector_Indices_With_Large_Residuals@commit` | 没有大残差时从非 void 函数末尾落出，是 UB；但其唯一上游 Greedy 本身不可达。 | T31 |
| 四模型 `eval_vector_interpolant_at_point@commit`（Vector Field 除外） | Single/Lajaunie/Stratigraphic/Continuous 创建 clone 后误用共享 kernel；顺序值可计算，并发调用不安全。 | T29 |
| `stratigraphic_surfaces.cpp::_get_increment_pairs@commit` | 同一对象二次 `ComputeInterpolant` 会保留旧 pair 元素并追加；新约束不一定进入实际 `[0]/[1]` pair。 | T26 |
| `continuous_property.cpp::get_equality_values@commit` | 添加被模型计数忽略的 planar/tangent 后仍按实际容器写 RHS，可越过按 interface 数量分配的 VectorXd。 | T27 |
| `modelling_input.cpp::spatial_metrics@commit` | 空约束公共查询返回 `resolution=0` 和极值 bounds，而非 `problem_computing_spatial_parameters`。 | T30 |

## 声明、定义与调用集合核对结论

- 30 个 T01 范围文件由首节四个精确的 `path-pattern@commit` 集合覆盖，并按表中
  符号证据再次核对；`geo_builder` 只作为整体排除边界。
- 头文件声明和 `.cpp`/内联定义双向扫描发现的唯一“非 pure virtual 且预期
  out-of-line、有声明无定义”函数是
  `Math_methods::_rot`；它也没有调用点，分类 `U`。
- 计划提示中的 `Point::operator==` 和注释中的 `outc`/GMP helpers 不在有效声明集，未伪造成能力。
- 具体核方法按上述四个笛卡尔积展开核对：完整定义 221 个，`R/AR` 哨兵 30 个；没有缺失
  `MQ3`、`WendlandC2` 或 `MaternC4` 定义。
- 五模型各 12 个虚槽逐一核对；所有具体覆写都有定义或内联体。显式 placeholder 恰为
  Stratigraphic Greedy 三项、Vector Greedy/reconstruction 四项、Continuous minimal/reconstruction
  两项，再加 R/AR 30 个导数和 coplanarity 一项。
- `TODO`、`TO implement`、空体、恒定返回、`throw -666`、throw 点和全仓调用点分别搜索；
  未把“测试未覆盖”当成不可达，`U` 均由零活动调用边证明。
- 每个表行只有一个分类码和一个关闭任务；展开集合无重叠。T02 只分类，没有实现、修改
  数学结果、创建新任务或宣称 oracle/parity/performance 通过。
