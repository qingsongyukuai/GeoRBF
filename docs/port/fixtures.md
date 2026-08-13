# Surfe golden fixture 与差分判定协议

## 身份、范围和非目标

golden fixture 唯一绑定冻结源码
`https://github.com/MichaelHillier/surfe.git@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
和 `georbf-surfe-oracle` protocol version 1。机器可读格式为
`tests/fixtures/schema/golden-fixture-v1.schema.json`，独立 Rust 读取、规范序列化和数值比较
基础位于 `tests/common/parity/mod.rs`。

T04 只固定格式、覆盖集、容差和判定流程；`tests/fixtures/schema/*.json` 是合成的 schema
正反例，不是 Surfe 数值输出，也不是正式 golden。正式 oracle 响应只能在相应迁移任务
完成实现和审阅后进入 fixture；T32 才关闭全局覆盖矩阵。本协议不实现核、矩阵、求解器或
模型算法，正常 Rust 测试不发现或启动 C++ oracle。

## 冻结源码依据

本协议的字段和分层判定对照了以下冻结实现；所有路径的提交均为
`290dbe0ab344f4258a4935f05cad0f153f0f69a4`：

- `surfe_lib/basis.{h,cpp}@commit`：`RBFKernel::{basis,dx_p1,dx_p2,dy_p1,
  dy_p2,dz_p1,dz_p2,dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz}`、
  `get_global_anisotropy`、`scaled_radius`、全部 `basis_*` 组合和
  `Modified_Kernel::basis_*`。核值、一阶导数、混合 Hessian、binary32 anisotropy 和
  Modified Kernel 必须属于不同数值层。
- `surfe_lib/matrix_solver.{h,cpp}@commit`：`Linear_LU_decomposition::solve`、
  `Quadratic_Predictor_Corrector::solve` 和
  `Quadratic_Predictor_Corrector_LOQO::solve`。分支、是否尝试、成功类别和迭代编号是精确
  字段；残差、可行性、目标、互补性和迭代浮点量分别按指定数值层比较。
- `surfe_lib/single_surface.{h,cpp}@commit`、`lajaunie.{h,cpp}@commit`、
  `stratigraphic_surfaces.{h,cpp}@commit`、`continuous_property.{h,cpp}@commit` 和
  `vector_field.{h,cpp}@commit`：各模型的 `get_interpolation_matrix`、
  `get_equality_values`、可达的 `get_inequality_matrix/get_inequality_values`、
  `setup_system_solver`、`eval_scalar_interpolant_at_point` 和
  `eval_vector_interpolant_at_point`。这些调用点固定矩阵/RHS、solver 分支、scalar 和
  gradient 的证据边界。

此处 `@commit` 是本节完整提交号的缩写；任务日志记录完整 `path@commit`。

## Fixture v1 envelope

每个 fixture 是一个 UTF-8 JSON 文件，恰好包含一个对象和结尾 LF。顶层字段及顺序固定为：

```text
schema, schema_version, fixture_id, source, generation, dataset,
comparison, request, expected
```

对象字段按 schema、oracle manifest 和本文件声明的顺序写出；数组保持输入或冻结源码
顺序。禁止时间戳、绝对路径、随机 id、地址、主机名以及未声明字段。`fixture_id` 和
`dataset.id` 是稳定的 ASCII 路径式 id，不包含平台或运行日期。

### 必需元数据

- `schema` 固定为 `georbf-surfe-golden`，`schema_version` 固定为 `1`。
- `source.repository`、`source.commit`、`source.oracle_protocol` 和 version 必须逐字匹配
  冻结身份。任何 commit 漂移都是 schema failure，不能按数值容差处理。
- `generation.command` 是 argv 数组，不是依赖 shell 解析的单个字符串。正式 fixture 必须
  记录实际 adapter 的显式生成命令；`schema-test-only` 仅可出现在 schema 合成正反例。
- `generation.environment` 固定记录 `OMP_NUM_THREADS=1`、`LC_ALL=C`、`TZ=UTC`。
- `generation.platform` 必须记录 OS、architecture、C++ compiler 完整版本、libc/运行时和
  endianness；不得记录绝对 reference 路径。
- `generation.precision` 固定为默认 `binary64`、anisotropy 中间量 `binary32`、
  `max_digits10` 和 `row_major`。平台差异依靠分层数值判定，不允许删除平台信息。
- `dataset.coverage` 使用 schema 的封闭枚举；`request_line_sha256` 和
  `response_line_sha256` 分别是规范 request/response JSON 行（包含末尾 LF）的 SHA-256。
- `request` 和 `expected` 是 oracle v1 request/response 对象。两者的 `request_id` 和
  `operation` 必须相同；expected source commit 必须再次精确匹配。

schema 表达 JSON 的结构约束；Rust validator 还拒绝重复键、错误字段顺序、不可表示为
有限 binary64 的普通数字、request/response 身份不一致、重复数值路径和非唯一数值分类。
两者都通过才是有效 fixture。

## 数值编码

有限量写作 JSON number，使用 locale-independent binary64 `max_digits10`。`-0.0` 必须
保留；fixture 再生成要求其字节精确，差分比较也要求零的符号一致。整数型离散量不得写成
浮点近似。

NaN 和无穷不得写作裸 `NaN`、`Infinity`、字符串或 JSON overflow；只允许以下单字段对象，
且比较时种类必须精确相同：

```json
{"number_kind":"nan"}
{"number_kind":"positive_infinity"}
{"number_kind":"negative_infinity"}
```

正常有效输入期望有限量。tagged non-finite 主要记录冻结错误/退化探针的可观察结果；它不把
失败转成成功。anisotropy 的 `storage: "f32"` 是精确字段，数值必须先在冻结 binary32
表达式中完成再扩展为 JSON number。

矩阵固定为 `{rows,cols,order,data}`，`order` 精确为 `row_major`，`data.length` 必须为
`rows * cols`。向量固定为 `{length,data}`，长度必须精确相等。形状、顺序、label、partition
和 polynomial term 都不是可容差字段。

## 默认 exact 与路径分类

每个 leaf 默认按 JSON 类型、对象字段、数组长度/顺序和值精确比较。只有
`comparison.numeric_rules` 中恰好匹配一个 `/` 开头路径模式的 leaf 才使用数值容差；`*`
只匹配单个 path segment。零个匹配即 exact，两个或更多匹配是 fixture 错误，不得选择较宽
规则。禁止 `all_numeric` 或其他单一全局容差。

必须精确的离散字段包括：

- schema/source/protocol/request identity、operation、coverage 和平台/精度标记；
- 核/模型/solver 名称和参数选择，约束接受或拒绝及 success/error 类别；
- 约束顺序和去重、level/group/reference/increment/stratigraphic pair 及源索引；
- 自由度计数、矩阵维数、row/column label、partition、row-major 顺序和 polynomial 项序；
- equality/inequality 符号分支、solver branch、`attempted`、`success`、active-set 离散选择；
- 迭代编号、迭代总数、停止原因和 error `stage/category/upstream_type/causes`。

`error.message` 只作诊断；它必须存在并在同一次生成的字节稳定性检查中稳定，但不参与
GeoRBF 行为 pass/fail。忽略列表只能显式写
`/expected/error/message` 或 `/expected/result/solve/weights/*`。不得把其他 mismatch 改成
diagnostic。

## 分层容差

有限数值使用包含边界的混合判定：

```text
delta = abs(actual - expected)
scale = max(abs(actual), abs(expected))
pass  = delta <= absolute + relative * scale
```

两值必须先为有限 binary64；tagged number 按上一节精确比较。下表是 fixture schema v1 的
唯一阈值，Rust `ToleranceClass::tolerance` 保存相同常量：

| class | absolute | relative | 用途 |
|---|---:|---:|---|
| `kernel_value` | `1e-12` | `1e-11` | 普通 isotropic/anisotropic 核值 |
| `first_derivative` | `1e-11` | `1e-10` | 对 point 1/2 的六个一阶导数 |
| `mixed_hessian` | `1e-10` | `1e-9` | `xx,xy,xz,yx,yy,yz,zx,zy,zz` |
| `anisotropy_f32` | `2e-6` | `2e-5` | transform、eigen-derived plunge、scaled radius 的 f32 路径 |
| `modified_kernel` | `2e-10` | `2e-9` | Lagrangian 修正后的全部泛函组合 |
| `matrix_rhs` | `1e-11` | `1e-10` | 完整 matrix、equality/inequality/range RHS |
| `solver_residual` | `1e-10` | `1e-8` | LU/QP residual、relative residual、complementarity |
| `solver_feasibility` | `1e-10` | `1e-8` | equality residual 与 inequality/bound violation |
| `prediction_scalar` | `1e-9` | `1e-8` | scalar field 和更新后的 iso-value |
| `prediction_gradient` | `1e-8` | `1e-7` | gradient field 三分量 |
| `solver_objective` | `1e-8` | `1e-7` | QP/LOQO objective |
| `iteration_numeric` | `1e-9` | `1e-7` | 逐迭代 residual、step、gap 等浮点证据 |

迭代 index/count 使用 exact，不属于数值 class。阈值不把非有限量变为可接受，也不替代
解析恒等式或有限差分门槛。T12–T14 的核/导数仍必须同时通过 C++ golden、解析性质和有限
差分三角验证。

## 权重、病态和非唯一系统

权重始终保存在 fixture 以便诊断和重建追踪，但不是独立的 release-blocking 等价条件。
`weight_policy` 只能是：

- `diagnostic_only`：不以权重 mismatch 判失败，适用于身份、错误或无求解用例；
- `residual_feasibility_predictions`：求解用例必须以实际尝试、有限性、残差、约束可行性、
  objective（QP/LOQO）及固定 witness points 的 scalar/gradient 共同判断。

`comparison.acceptance` 机器可读地保存 `required_finite`、`residual_l2_max`、
`relative_residual_max`、`equality_residual_linf_max`、`inequality_violation_linf_max` 和
`prediction_witnesses_required`。无求解用例的四个 ceiling 必须为 `null`；求解成功用例不得
用 `null` 省略实际适用的 residual/feasibility 门槛。

病态或非唯一系统不得因 condition number 提前拒绝。fixture 必须精确比较 solver branch 和
`attempted: true`，然后要求：所有要求有限的证据有限；residual 和 feasibility 分别满足
其数值层以及该数据集记录的绝对可行性 ceiling；QP/LOQO 的 objective 和 complementarity
满足对应层；所有 witness prediction 满足 scalar/gradient 层。只有这些同时满足才能通过。
权重相同不能弥补 residual、feasibility 或预测失败，权重不同也不能单独制造失败。

每个求解 fixture 必须记录与数据尺度相关的 residual/feasibility ceiling；该 ceiling 是
额外的单边门槛，不得比上述 golden 差分层更宽来掩盖 mismatch。不可行、奇异或非有限
用例按 error stage/category 精确比较，不能用宽容差转成成功。

## 确定性覆盖数据集

以下是固定 case family 清单。坐标、方向、shape 和 RHS 只能使用已列的十进制常量或由其
确定性笛卡尔积产生；数组按这里的书写顺序输出。此清单定义将来 fixture id 和必须证据，
不包含 oracle 输出。

### Kernel families

- `kernel/isotropic/<kernel>/separated`：九核 `cubic,gaussian,mq,mq3,tps,imq,r,
  wendland_c2,matern_c4`；point 1 `[-1.25,0.5,2.0]`、point 2
  `[0.75,-0.25,1.0]`、shape `0.7`，覆盖 basis、point 1/2 first derivatives 和完整 mixed
  Hessian。
- `kernel/isotropic/<kernel>/zero-near-support`：同点、`2^-20` 轴向间隔、负坐标、多尺度
  `1e-3/1/1e3`；compact/support 核另测 `radius = support * (1-2^-20), support,
  support * (1+2^-20)`。分支结果和 tagged non-finite 精确。
- `kernel/functionals/<kernel>/directions`：direction 1 `[0.3,-0.4,0.5]`、direction 2
  `[-0.2,0.7,0.1]`，覆盖 Value/Planar/Tangent 的全部 oracle `basis_*` 组合和参数交换。
- `kernel/anisotropy/identity-oblique-degenerate`：单位轴 normals、oblique normals
  `[0.36,-0.48,0.8]` 的固定排列、两个 planar 下限和 eigenvalue `0.0001` 截断邻域；保存
  `storage=f32`、transform、global plunge、scaled radius、核与导数。
- `kernel/modified/<kernel>/all-functionals`：unisolvent candidates 依次为
  `[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1]`，覆盖选中源索引、16 个 Lagrangian 系数和
  Value/Planar/Tangent 全组合；另有 insufficient/coplanar error case。

### Matrix and model families

模型 family 共用下列有序、精确十进制数据；某模型不接受的约束类别从请求中省略，不以零值
占位：

- `single-base`：interface `[x,y,z,level]` 依次为
  `[0,0,0,0],[1,0,0,0],[0,1,0,0],[1,1,0,0]`；planar
  `[x,y,z,nx,ny,nz] = [0.5,0.5,0,0,0,1]`；tangent
  `[x,y,z,tx,ty,tz] = [0.25,0.25,0,1,0,0]`；inequality
  `[x,y,z,sign]` 依次为 `[0.5,0.5,0.75,1]`、`[0.5,0.5,-0.75,-1]`；evaluation points
  依次为 `[0.25,0.75,0.5]`、`[0.75,0.25,-0.5]`、`[0.5,0.5,0]`。
- `multilevel-base`：interface 依次为
  `[0,0,-1,-1],[1,0,-1,-1],[0,0,0,0],[0,1,0,0],[0,0,2,2],
  [1,1,2,2]`；planar 为 `[0.5,0.5,0,0,0,1]`；evaluation points 依次为
  `[0.25,0.25,-0.5]`、`[0.5,0.5,1]`、`[0.75,0.75,2.5]`。输入顺序故意不按 level
  二次重排；只有冻结源码排序可以改变响应顺序。
- `continuous-base`：interface 依次为
  `[0,0,0,-1],[1,0,0,0.5],[0,1,0,1.25],[0,0,1,2]`，planar 为
  `[0.25,0.25,0.25,0.2,-0.4,0.8]`，evaluation points 使用
  `[0.2,0.3,0.1]`、`[0.6,0.1,0.2]`、`[-0.2,0.4,0.5]`。
- `vector-base`：planar 依次为 `[0,0,0,1,0,0]`、`[1,0,0,0,1,0]`、
  `[0,1,1,0.36,-0.48,0.8]`，evaluation points 使用 `[0.25,0.25,0.25]`、
  `[0.75,-0.25,0.5]`、`[-0.5,0.5,1.25]`。

- `model/single_surface/equality`、`inequality-active-inactive`、`restricted-range`：固定四个
  interface、一个 planar、一个 tangent、上下 inequality 和三个 evaluation points，全部取
  `single-base`；覆盖 normalized constraints、layout、完整 matrix/RHS、LU/QP/LOQO、scalar
  和 gradient。
- `model/lajaunie/multilevel-increments`：levels `-1,0,2`，每层一/多点，固定 reference 和
  increment pairs，输入取 `multilevel-base`；另含 restricted-range 与 iso-value update 证据。
- `model/stratigraphic/three-level-lithology`：三个精确 levels、层内 difference、边界层和
  层间/层外 inequality，输入取 `multilevel-base` 并按顺序追加 inequality
  `[0.5,0.5,-1.5,-1]`、`[0.5,0.5,0.75,1]`、`[0.5,0.5,2.5,1]`；覆盖上下邻层 pair、
  最小层间约束、matrix/RHS 和最终层位关系。
- `model/continuous_property/reachable`：只请求冻结公开 API 可达的 interface/planar 数据、
  matrix/RHS、scalar/gradient，输入取 `continuous-base`；TODO 或不可达分支使用 error case，
  不补写能力。
- `model/vector_field/planar-hessian`：三个非共线 planar normals 和三个 evaluation points；
  输入取 `vector-base`，覆盖 component layout、Hessian matrix、normal RHS、potential scalar
  和 gradient。

上述 model family 各自至少包含无 polynomial、0/1/2 阶中源码可达的代表项、smoothing 关闭
和一个非零 smoothing 用例；所有新增组合只能归入这些固定 family，不改变任务序列。

### Solver and error families

- `solver/lu/well-conditioned`、`ill-conditioned-attempted`、`singular`、`non-finite`：显式
  row-major matrix/RHS；well-conditioned 使用 `[[4,1],[1,3]]/[1,2]`，ill-conditioned 使用
  3×3 Hilbert matrix（元素按 binary64 `1/(row+col+1)` 生成）且 RHS 精确由全一向量乘出，
  singular 使用 `[[1,2],[2,4]]/[3,6]`，non-finite 以 tagged value 替换第一个 matrix 元素；
  保存 attempted、pivot/branch、weights、residual 和 failure 类别。
- `solver/qp/equality-inequality`、`active-boundary`、`inactive-boundary`、`infeasible`：固定
  Hessian `[[2,0],[0,2]]`、equality row `[[1,1]]`/RHS `[1]`、inequality rows
  `[[1,0],[0,1]]`；lower RHS 依次用 `[0,0]`、`[0.75,0]`、`[-1,-1]`，infeasible 追加
  `[-1,0] >= 0` 与 `[1,0] >= 1`；保存 objective、feasibility、complementarity 和逐迭代证据。
- `solver/loqo/single-double-tight-bound`、`infeasible`：固定 lower/range 语义，保存与普通
  QP 相同 Hessian，`A` 为 identity；`b/r` 依次使用 `[0,0]/[1,1]`、
  `[0.25,0.25]/[0.5,0.5]`、`[0.5,0.5]/[0,0]`，不可行 case 使用负 range；保存与普通 QP
  同级的证据和 restricted-range branch。
- `error/<stage>/<category>`：覆盖 request、configuration、constraint ingest、preprocess、
  basis、assembly、solve、reconstruction、evaluation 和 oracle safety；stage/category/causes
  精确，message 仅诊断。

每个正式 fixture 必须把所需 evidence 逐项列入 `dataset.coverage`；不能靠 case 名声称覆盖。
T32 的覆盖报告必须证明上述每个 family 已生成、审阅和执行，没有 skip。

## 生成、重复性和审阅

正式生成必须先按既定优先级解析 reference，确认 HEAD 精确且工作树干净，然后确认 oracle
identity。推荐流程如下；实际 argv 和平台版本写入 fixture：

```sh
git -C "$resolved_reference" rev-parse HEAD
git -C "$resolved_reference" status --short
OMP_NUM_THREADS=1 LC_ALL=C TZ=UTC .cache/surfe-oracle/surfe-oracle \
  < .cache/surfe-oracle/request.jsonl \
  > .cache/surfe-oracle/response-1.jsonl
OMP_NUM_THREADS=1 LC_ALL=C TZ=UTC .cache/surfe-oracle/surfe-oracle \
  < .cache/surfe-oracle/request.jsonl \
  > .cache/surfe-oracle/response-2.jsonl
cmp .cache/surfe-oracle/response-1.jsonl \
    .cache/surfe-oracle/response-2.jsonl
```

两次响应必须逐字节相同；随后对 request/response 规范行（含 LF）计算 SHA-256，组装 fixture，
解析后按 Rust canonical serializer 再写出两次并要求字节相同。任何 compiler/platform、请求、
oracle adapter 或 source commit 变化都必须重新生成受影响 fixture，不能手改 expected 数字。

审阅者必须同时检查：request 的地质/数学意图；source 和平台身份；coverage 是否与 response
证据一致；每个 numeric path 恰有一个正确 class；所有其他 leaf 是否 exact；matrix/vector
长度与 label；solver 是否实际尝试；非唯一权重是否使用 residual/feasibility/prediction；两次
输出和两次 serialization 是否字节稳定。正式 fixture 更新必须与产生它的实现/parity 任务
一起提交，并在日志记录生成命令、hash、比较结果和审阅人/方式；未经审阅的输出不得提交。

## T04 验证入口

当前仓库尚无 Cargo 工程，因此本任务运行：

```sh
rustc --edition=2021 --test tests/common/parity/mod.rs -o /tmp/georbf-parity-tests
/tmp/georbf-parity-tests
```

该测试覆盖正反例识别、规范序列化往返、容差边界、signed zero、tagged non-finite、路径唯一
分类和分层阈值。JSON Schema 还要由 Draft 2020-12 validator 对一个正例和两个反例验证。
一旦 Cargo 工程在 T05 建立，此模块应进入普通 `cargo test`，但不得增加 oracle 运行依赖。
