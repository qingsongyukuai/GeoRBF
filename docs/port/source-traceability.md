# GeoRBF 逐模块来源追踪

## 冻结身份与使用方式

所有行为来源均固定到 Surfe 提交
`290dbe0ab344f4258a4935f05cad0f153f0f69a4`。下表逐一覆盖当前 43 个
生产 Rust 模块；“来源与符号”描述行为证据，不表示仓库包含 C++ 文件。更细的逐任务
阅读、oracle、差分与安全处理证据在 `JOURNAL.md`，符号覆盖基线在 `inventory.md`、
`symbol-classification.md` 和 `parity-report.md`。

## 生产模块映射

| GeoRBF 模块 | 冻结来源与主要符号 |
|---|---|
| `src/lib.rs` | `surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：安全公开入口与五模型集成边界。 |
| `src/parameters.rs` | `surfe_lib/modelling_parameters.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：枚举、默认、名称与 setter。 |
| `src/error.rs` | `surfe_lib/grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：23 个异常类别。 |
| `src/geometry.rs` | `surfe_lib/modelling_input.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：Point 字段与安全值边界。 |
| `src/constraints.rs` | `surfe_lib/modelling_input.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：Interface/Inequality/Planar/Tangent 与方向转换。 |
| `src/constraints/grouping.rs` | `surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：exact level、reference 与 interface grouping。 |
| `src/ordering.rs` | `surfe_lib/modelling_input.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`math_lib/math_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：排序、collocation、indexed sort。 |
| `src/spatial.rs` | `surfe_lib/modelling_input.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：距离、邻居、bounds、extremal 与 spatial metrics。 |
| `src/polynomial.rs` | `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：0/1/2 阶 polynomial 与 truncation。 |
| `src/polynomial/lagrangian.rs` | `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：unisolvent 选择与 Lagrangian basis。 |
| `src/kernel/mod.rs` | `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：核类型集成边界。 |
| `src/kernel/derivatives.rs` | `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：一阶/混合二阶分量派发。 |
| `src/kernel/isotropic.rs` | `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：九种 isotropic 核、两点一阶导数与 mixed Hessian。 |
| `src/kernel/anisotropy.rs` | `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：global anisotropy、scaled radius 与六种 anisotropic 核。 |
| `src/kernel/modified.rs` | `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：Modified_Kernel 全组合。 |
| `src/functional.rs` | `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 及五模型 `get_interpolation_matrix`：Value/Derivative/Tangent/Difference。 |
| `src/layout.rs` | `surfe_lib/{single_surface,lajaunie,stratigraphic_surfaces,continuous_property,vector_field}.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：五模型自由度顺序与分区。 |
| `src/assembly.rs` | `surfe_lib/{single_surface,lajaunie,stratigraphic_surfaces,continuous_property,vector_field}.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：matrix/RHS/polynomial/smoothing。 |
| `src/solver/mod.rs` | `surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`math_lib/math_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：求解器集成边界。 |
| `src/solver/error.rs` | `surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：LU 失败阶段与外层类别。 |
| `src/solver/lu.rs` | `surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：Linear_LU_decomposition；`Eigen/src/LU/PartialPivLU.h@36b95962756c1fce8e29b1f8bc45967f30773c00` 仅作冻结求解路径证据。 |
| `src/solver/qp_predictor_corrector.rs` | `surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`math_lib/math_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：ordinary predictor-corrector QP。 |
| `src/solver/qp_loqo.rs` | `surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`math_lib/math_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：restricted-range/LOQO QP。 |
| `src/model/mod.rs` | `surfe_lib/modeling_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、五模型头文件：模型模块边界。 |
| `src/model/fitted.rs` | `surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：ComputeInterpolant、getter 与单点/批量 field。 |
| `src/model/reconstruct.rs` | `surfe_lib/{single_surface,lajaunie,stratigraphic_surfaces}.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：三个 convert_modified_kernel_to_rbf_kernel body。 |
| `src/model/single_surface/layout.rs` | `surfe_lib/single_surface.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：Single layout 与分区。 |
| `src/model/single_surface/assembly.rs` | `surfe_lib/single_surface.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：Single equality/inequality/bounded values。 |
| `src/model/single_surface/mod.rs` | `surfe_lib/single_surface.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/modeling_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：linear/QP/restricted fit 与 field。 |
| `src/model/lajaunie/layout.rs` | `surfe_lib/lajaunie.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：reference/increment layout。 |
| `src/model/lajaunie/assembly.rs` | `surfe_lib/lajaunie.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：equality 与 restricted bounds。 |
| `src/model/lajaunie/mod.rs` | `surfe_lib/lajaunie.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/modeling_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：fit、iso update 与 field。 |
| `src/model/stratigraphic/layout.rs` | `surfe_lib/stratigraphic_surfaces.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：层序、lithology 和 same-level layout。 |
| `src/model/stratigraphic/assembly.rs` | `surfe_lib/stratigraphic_surfaces.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：ordinary/restricted matrix partition、RHS 与 bounds。 |
| `src/model/stratigraphic/mod.rs` | `surfe_lib/stratigraphic_surfaces.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：层序拟合、求解、conversion、iso 与 field。 |
| `src/model/continuous_property/layout.rs` | `surfe_lib/continuous_property.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：实际可达 interface layout。 |
| `src/model/continuous_property/assembly.rs` | `surfe_lib/continuous_property.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：实际可达 equality RHS。 |
| `src/model/continuous_property/mod.rs` | `surfe_lib/continuous_property.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：可达 fit/field 与 TODO/UB 边界。 |
| `src/model/vector_field/layout.rs` | `surfe_lib/vector_field.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：Planar 三分量 layout。 |
| `src/model/vector_field/assembly.rs` | `surfe_lib/vector_field.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：normal RHS。 |
| `src/model/vector_field/mod.rs` | `surfe_lib/vector_field.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：Hessian fit、potential 与 gradient。 |
| `src/builder.rs` | `surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：Add/Set/Get、参数与安全 fit lifecycle。 |
| `src/greedy.rs` | `surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/modeling_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、五模型 hook：公开零轮可达行为。 |

## 非生产证据与许可证

- `tests/fixtures/golden/global-parity-v1.json` 保存冻结 oracle transcript、hash 和 family
  索引，不包含 C++ 源码；生成和覆盖见 `parity-report.md`。
- `benches/` 对应冻结五模型预处理、装配、求解与批量评估热路径；公平性和 adapter
  身份见 `performance-report.md`。
- 外部 Eigen 路径只用于复现冻结 binary32 eigensolver/partial-pivot LU 的可观察次序，
  不被 vendored、链接或打包；GeoRBF Cargo 依赖为零。
- Rust 翻译和改编所需的 Surfe 2017 Government of Canada MIT notice 随根 `NOTICE`
  进入发布包；项目自身许可证为根 `LICENSE`。
