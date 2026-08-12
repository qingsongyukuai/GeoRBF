# 冻结 Surfe reference oracle 协议

## 身份、用途与边界

本协议唯一绑定
`https://github.com/MichaelHillier/surfe.git@290dbe0ab344f4258a4935f05cad0f153f0f69a4`。
oracle 用冻结 C++ 生成后续差分测试所需的可观察结果和中间证据，不是 GeoRBF 的生产
实现，也不定义新的数学语义。机器可读字段权威清单位于
`docs/port/oracle-manifest.json`。

reference 仍按 `SURFE_REFERENCE_DIR`、`../surfe`、`.cache/surfe-reference` 的唯一顺序
发现；oracle 源适配层、对象、库、可执行文件、请求和响应只能位于仓库外或被忽略的
`.cache/surfe-oracle`。正常 Cargo 构建、测试、文档测试和发布不得发现、构建或调用
oracle。不得提交 Surfe、Eigen、C++ adapter、对象、共享库、可执行文件或本任务的
数值响应。

adapter 只能做以下三类工作：把协议输入转换为冻结类型；调用冻结函数；只读提取并
序列化结果。不得更改公式、分支、排序、矩阵布局、solver 参数、停止规则或异常行为。
需要读取 `protected`/`private` 状态时，探针必须只存在于外部 adapter，且只能读取；
任何用于追踪局部迭代变量的插桩都必须在外部副本中保持表达式和执行顺序不变，并在
对应后续任务同时核对未插桩最终输出。T03 smoke 没有修改 reference 文件。

## T03 审阅的冻结源码与提取点

下列路径均带固定提交
`290dbe0ab344f4258a4935f05cad0f153f0f69a4`：

- `surfe_lib/surfe_api.{h,cpp}@commit`：约束 `Add*/Set*`、参数 setter、
  `ComputeInterpolant`、`EvaluateInterpolantAtPoint(s)`、
  `EvaluateVectorInterpolantAtPoint(s)`、约束 getter 和公开错误边界。
- `surfe_lib/basis.{h,cpp}@commit`：`RBFKernel::{basis,dx_p1,dx_p2,dy_p1,
  dy_p2,dz_p1,dz_p2,dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz}`、
  `get_global_anisotropy`、`scaled_radius`、全部 `basis_*` 泛函组合、
  `Lagrangian_Polynomial_Basis` 和 `Modified_Kernel::basis_*`。
- `surfe_lib/modeling_methods.{h,cpp}@commit`：约束清洗、interface level/group/reference、
  `setup_basis_functions`、`get_equality_matrix`、solver/kernel/weights 所有权和错误包装。
- `surfe_lib/matrix_solver.{h,cpp}@commit`：`Linear_LU_decomposition`、
  `Quadratic_Predictor_Corrector`、`Quadratic_Predictor_Corrector_LOQO` 的输入、分支、
  `solve`、验证和最终 `weights`。
- `surfe_lib/single_surface.{h,cpp}@commit`、`lajaunie.{h,cpp}@commit`、
  `stratigraphic_surfaces.{h,cpp}@commit`、`continuous_property.{h,cpp}@commit`、
  `vector_field.{h,cpp}@commit`：各自 `process_input_data`、`get_method_parameters`、
  `get_interpolation_matrix`、equality/inequality RHS、solver 设置和 scalar/vector 评估；
  Lajaunie/Stratigraphic 的 increment 与层序 pair 由外部只读探针记录。
- `surfe_lib/modelling_input.{h,cpp}@commit`：本任务构造请求和序列化清洗后约束所需的
  `Point/Interface/Planar/Tangent/Inequality/Constraints` 字段。
- `math_lib/math_methods.{h,cpp}@commit`：QP/LOQO 目标、步长、残差和迭代证据的实际
  计算来源；只在 `solver.run` 请求该证据时读取。
- `test/main.cpp@commit`：公开构建/求解调用顺序的使用证据；其中 `Geo_Builder`、VTK、
  Qt、网格、等值面写出和可视化不进入 oracle。

路径中的 `@commit` 是上段完整提交的缩写，仅为提高表格可读性；日志使用完整提交号。

## JSON Lines 传输

stdin 每行恰好一个 UTF-8 JSON request，stdout 对每个 request 恰好输出一个 JSON
response 和一个 LF。正常 stdout 不得混入 Surfe 进度、condition number 或调试文本；
诊断只写 stderr。T03 固定为一次进程处理一个请求，后续可以保持相同 envelope 扩展为
多行流，但不得改变一请求一响应关系。

request 顶层字段顺序固定为：

```text
protocol, protocol_version, request_id, source_commit, operation, input, evidence
```

response 顶层字段顺序固定为：

```text
protocol, protocol_version, request_id, source, operation, status, result, error
```

`protocol` 必须为 `georbf-surfe-oracle`，`protocol_version` 必须为 `1`，请求和响应的
提交必须逐字等于固定提交。未知顶层字段、未知 operation、提交不匹配、重复 request
id 或非对象 `input` 必须返回 request-stage error，不能猜测或忽略。

`status: "ok"` 时 `result` 必须是对象且 `error` 必须为 `null`；`status: "error"`
时正好相反。错误对象的 `category` 是后续 parity 比较字段，`message` 只是诊断文字。
嵌套异常按外到内放入 `causes`，不得用异常消息文本代替稳定类别。

## 数值、矩阵和确定性

- 普通数值用 locale-independent binary64 JSON number，输出精度至少为 C++
  `max_digits10`。负零写为 `-0.0`；NaN 和正负无穷使用 manifest 的 tagged object，绝不
  输出无效 JSON token。
- Surfe 的 binary32 anisotropy 中间量必须先在冻结表达式中完成 `f32` 运算，再提升为
  JSON number，并在容器中记录 `storage: "f32"`；不得在 oracle 中改成全 `f64`。
- 矩阵总是 `{rows, cols, order: "row_major", data}`，即使 Eigen 内部默认列主序；
  `data.length` 必须精确等于 `rows * cols`。向量总是 `{length, data}`。
- 约束和索引数组保持冻结源码顺序；只有冻结实现本身排序时才排序。对象字段按 manifest
  的 `ordered_fields` 输出，不包含绝对路径、时间戳、随机 id 或地址。
- 每次运行固定 `OMP_NUM_THREADS=1 LC_ALL=C TZ=UTC`。相同 executable、相同请求、相同
  环境必须产生逐字节相同 stdout；不同编译器/平台的数值比较阈值由 T04 定义，T03 不
  预先设置宽松容差。

## Operation 与证据字段

`identity` 返回固定仓库/提交、编译信息和 adapter 支持的 operation/evidence，不运行
算法。其结果用于拒绝连接到漂移的 reference。

`kernel.evaluate` 接受核名、shape、两点、可选方向与 anisotropy/interface groups。
它可请求：核值；分别对 point 1/2 的六个一阶导数；按 `xx,xy,xz,yx,yy,yz,zx,zy,zz`
排列的混合 Hessian；全局 anisotropy 的 binary32 transform/plunge/scaled radius；
Value、Planar、Tangent 的全部 `basis_*` 组合；Modified Kernel 的 unisolvent 原始索引、
Lagrangian 系数和全部组合。`R/AR` 的冻结 `-666` 哨兵按实际结果编码，不能由 adapter
替换为“正确”导数。

`model.run` 接受五种模型之一、参数、四类约束和有序评估点。可请求：清洗后的完整
约束；精确 interface level/group/reference/increment/stratigraphic pair；带行列标签、
分区和 polynomial 项的 layout；完整 interpolation/equality/inequality matrix；所有
RHS/range；solver 证据；按输入顺序的 scalar/gradient。矩阵和 RHS 必须通过各模型
公开虚函数提取，不能在 adapter 重新实现装配公式。

`solver.run` 直接接收显式矩阵/RHS，选择 partial-pivot LU、普通 predictor-corrector QP
或 LOQO restricted-range 分支。必须记录“是否实际尝试”、最终有限性、权重、残差、
可行性、目标、互补性和可获得的迭代证据；不得根据 condition number 提前拒绝。

`error.probe` 从指定公开入口和状态执行预期失败。错误阶段至少区分 request、配置、
约束摄取、预处理、basis、装配、求解、重建和评估。C++ 未定义行为不执行；这类用例
返回稳定的 oracle safety error，并由对应后续兼容任务决定 Rust 类别。

## 外部构建与运行

从 GeoRBF 根目录解析 `resolved_reference` 后，先执行：

```sh
git -C "$resolved_reference" rev-parse HEAD
git -C "$resolved_reference" status --short
git -C "$resolved_reference" submodule update --init eigen-git-mirror
git -C "$resolved_reference" submodule status eigen-git-mirror
```

HEAD 必须为固定提交，工作树必须干净，Eigen gitlink 必须为
`36b95962756c1fce8e29b1f8bc45967f30773c00`。不要初始化或构建 pybind11、Qt、VTK 或
`geo_builder`。在 `.cache/surfe-oracle` 放置经审阅、实现本 manifest 的外部 adapter；
编译时只包含冻结的 `surfe_lib`、`math_lib` 和该 Eigen gitlink，并使用 release/
`NDEBUG` 关闭调试 stdout。T03 实测使用 `g++ -std=c++11 -O2 -DNDEBUG -fopenmp`，编译
`math_methods.cpp`、九个非可视化核心 `.cpp`、`surfe_api.cpp` 与外部 adapter；adapter
文件和产物均被忽略。

调用和确定性检查形式固定为：

```sh
OMP_NUM_THREADS=1 LC_ALL=C TZ=UTC .cache/surfe-oracle/surfe-oracle \
  < .cache/surfe-oracle/request.jsonl \
  > .cache/surfe-oracle/response-1.jsonl
OMP_NUM_THREADS=1 LC_ALL=C TZ=UTC .cache/surfe-oracle/surfe-oracle \
  < .cache/surfe-oracle/request.jsonl \
  > .cache/surfe-oracle/response-2.jsonl
cmp .cache/surfe-oracle/response-1.jsonl \
    .cache/surfe-oracle/response-2.jsonl
```

之后用 `jq` 校验 envelope、提交、operation、status、矩阵/向量长度和 operation-specific
必需字段，并可把完整 response 行通过环境变量交给
`tests/common/oracle_protocol.rs` 的外部 smoke 校验。响应只能留在 ignored cache；正式
golden fixture 的生成、容差和审阅从 T04 开始。

## T03 实测 smoke

本任务在 `.cache/surfe-oracle` 构建外部 adapter，运行
`model.run/single_surface_linear_smoke`：四个 interface 点、一个 planar 点、Cubic 核、
一阶 polynomial，输出 11×11 matrix、11 项 RHS/LU weights、残差和一个 scalar/gradient
预测；公开 API 结果同时与直接模型证据路径交叉核对。相同请求连续两次的响应
SHA-256 均为
`b15aaba2405d51278d61d1d4b0d5e57c0c4fcb2ae1adc1cecabfbd067321cc7c`，`cmp` 通过；
响应内 commit 精确匹配、所有数值有限、LU `attempted/success` 为真、`residual_l2 = 0`。
这些是 oracle 可运行和协议确定性的证据，不是正式 golden fixture，也不代表 T32 全局
parity 或 T33 性能通过。
