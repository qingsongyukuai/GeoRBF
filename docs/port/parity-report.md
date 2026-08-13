# T32 全局差分 parity 报告

## 结论

冻结 Surfe 非可视化核心在 T00–T31 确认纳入且实际可达的离散行为和数值行为，已由正式
fixture、冻结 C++ transcript、Rust golden tests、解析恒等式与有限差分共同关闭。T32 的
release-blocking parity 结果为 **PASS**：没有未解释 mismatch、没有 `#[ignore]`、没有未运行的
固定 family，也没有通过放宽 T04 容差、删除 fixture 或改变求解定义制造通过。

本结论只覆盖行为 parity。T33 才执行同机性能比较；本报告不宣称性能达到或超过 Surfe。

## 冻结身份与源码复核

- 源仓库：`https://github.com/MichaelHillier/surfe.git`。
- 固定提交：`290dbe0ab344f4258a4935f05cad0f153f0f69a4`。
- reference：按规定优先级解析到 `.cache/surfe-reference`；T32 开始和结束均以
  `git rev-parse HEAD` 得到固定提交，`git status --short` 为空。
- Eigen gitlink：`36b95962756c1fce8e29b1f8bc45967f30773c00`，只用于仓库外冻结 oracle；
  Cargo 构建、测试和发布不依赖它。

T32 重新枚举并校验了 28 个冻结 `.h/.cpp` 文件，定义/调用/缺陷审计得到 631 条证据，审计
文件 SHA-256 为 `1f0f8ec66e5a24071035f7ed5ab770071943743edf8fd9d7cfedb4449a8ccec4`。
复核的准确来源与主要符号如下；每个路径均位于上述固定提交：

| 来源 | T32 复核符号/边界 |
|---|---|
| `surfe_lib/modelling_parameters.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 全部参数、枚举、默认与 23 个错误类别 |
| `surfe_lib/modelling_input.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 约束、方向、清洗、排序、空间与 residual selector |
| `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | polynomial、Lagrangian、九核、两点一阶导数、mixed Hessian、anisotropy、Modified Kernel 和全部 `basis_*` |
| `math_lib/math_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | indexed sort、QP/LOQO、步长 helpers 与 `partialPivLu` 调用 |
| `surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | LU、ordinary predictor-corrector、restricted LOQO 的验证/尝试/停止/失败 |
| `surfe_lib/modeling_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | grouping、basis setup、equality slicing、工厂、iso update、Greedy 调用集 |
| `surfe_lib/single_surface.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | linear、ordinary QP、restricted、matrix/RHS、field 与 conversion body |
| `surfe_lib/lajaunie.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | reference、same-level increment、iso update、restricted 与 field |
| `surfe_lib/stratigraphic_surfaces.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 层序、上下邻层、lithology pairs、QP/LOQO、conversion 与 field |
| `surfe_lib/continuous_property.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 公开可达 interface value 系统、不可达/TODO 与越界缺陷边界 |
| `surfe_lib/vector_field.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | Planar Hessian、normal RHS、势函数与梯度 |
| `surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | 全部 Add/Set/Get、`ComputeInterpolant`、单点/批量 scalar/gradient |
| `surfe_lib/{debug,surfe_lib_module}.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`math_lib/math_lib_module.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` | debug/visualization 排除和 DLL/CMake 导出边界 |

定义/调用集合与 `docs/port/inventory.md`、`symbol-classification.md` 的纳入、不可达、缺陷和
排除分类一致；T32 没有发现新的未归属可达符号。

## 正式 fixture 与可重复生成

正式 fixture 为 `tests/fixtures/golden/global-parity-v1.json`，遵循
`georbf-surfe-golden` v1 和 `georbf-surfe-oracle` v1。它不是只含 checksum 的空索引：
`expected.result.capabilities.probe_catalog` 保存 T06–T31 共 26 个冻结 C++ probe 的完整 stdout/
stderr transcript、调用参数、退出码、source/binary/transcript SHA-256；
`family_coverage` 保存 74 个固定/补充 family 到 Rust release-blocking test 的唯一交叉索引。

生成环境固定为 `OMP_NUM_THREADS=1`、`LC_ALL=C`、`TZ=UTC`；实测平台为 Linux x86_64、
g++ 13.3.0、glibc 2.39。被忽略生成器
`.cache/surfe-oracle/t32_generate_fixture.py` SHA-256 为
`d33d54e4b3ed2fd4d71020059268f1bfc4eb5cdc92939c304fd9c2fe55bdcb53`。生成器先按规定优先级
验证 reference HEAD 和 clean status，再执行每个 probe invocation；每个 invocation
在一次生成内运行两遍并要求字节相同；整个 fixture 又独立生成两遍并由 `cmp` 确认相同。
最终 fixture 为 108,355 bytes，SHA-256 为
`6c26cb6bdd1f2bdba00f06c9dcdb27415cc45ecb1fb4d4e1569cc1fc3a754c28`。

正常 `cargo test` 不发现、不启动 C++ oracle。`tests/parity/main.rs` 只读取已审阅 fixture，验证：

- schema、固定提交、字段顺序和 canonical byte round-trip；
- request/response 行 SHA-256；
- 每个保存 transcript 的 SHA-256；
- T06–T31 无缺失/重复 probe；
- 74 个 family 无缺失/重复，且每项指向实际存在的 Rust test function；
- summary 中 ignored tests 和 unexplained mismatches 均为零。

## 覆盖矩阵

机器可读逐 family 明细位于正式 fixture；以下是审阅汇总。

| 固定范围 | family 数 | C++ 证据 | Rust/三角证据 | 结果 |
|---|---:|---|---|---|
| 九种 isotropic 核：separated、zero/near/support | 18 | T12 transcript，value/六个一阶/九个 Hessian 槽及分支 | exact bits、交换/符号/对称恒等式、两点有限差分 | PASS |
| Value/Derivative/Tangent 与 Modified Kernel | 18 | T14/T15 transcript，25 个 primitive、difference 和全部 modified 组合 | 泛函线性/交换、Lagrangian 消去、两层有限差分；九核基面由 T12 exact golden 约束 | PASS |
| binary32 全局各向异性 | 1 | T13 eigen/transform/plunge/kernel transcript | f32 bit evidence、截断分支、各向同性极限、有限差分 | PASS |
| Single equality / QP / restricted | 3 | T22–T24 完整 matrix/RHS/solver/field transcript | layout exact、residual/feasibility/objective、scalar/gradient witness | PASS |
| Lajaunie / Stratigraphic / Continuous / Vector | 4 | T25–T28 模型 transcript | reference/pair/layer exact、matrix/RHS、solver 后验、fields | PASS |
| LU / ordinary QP / LOQO | 10 | T18–T20 正常、病态、奇异、非有限、活跃、不可行和 tight transcript | attempted/branch/iteration exact；residual、feasibility、objective、complementarity | PASS |
| request 至 oracle-safety 的错误阶段 | 10 | T27/T30 公开错误与安全探针 | 类型化类别 exact；message 仅诊断；UB/terminate 不复制 | PASS |
| geometry、ordering、spatial、polynomial、layout、assembly、reconstruction、API、Greedy | 10 | T07–T11、T16–T17、T21、T29–T31 transcript | 离散 exact、matrix bit hash、public field、零轮 Greedy | PASS |

九核 family 中 `R` 的 value 与错误 sentinel 都作为离散行为验收；不存在的导数没有被伪造成
数值能力。anisotropic 工厂只接受冻结源码具有类实现的六种核。病态或非唯一系统先实际求解，
权重只作诊断；通过条件始终是有限性、残差、约束可行性、objective/complementarity 和固定
witness prediction，未使用 condition-number pre-gate。

## 已解释差异与缺陷兼容

T32 复跑没有发现新 mismatch。既有差异均已由兼容测试和 `docs/port/compatibility.md` 关闭：

- TPS/Wendland 非零第四坐标 Hessian、Modified Kernel 固定 Gram 对角和 Hessian 乘法次序；
- Eigen/Rust 病态或非唯一权重不同但 residual/feasibility/prediction 等价；
- ordinary QP/LOQO 的冻结错误成功、非终止或不合格 terminal candidate，在保存 source
  candidate/trace 后由 Rust 类型化后验安全拒绝；
- `SetTangentConstraints` 错写 Planar、Continuous Property planar/tangent 越界、未初始化状态、
  nested exception terminate 等 C++ UB/缺陷，不复制到 Rust；
- conversion bodies 和 Greedy hooks 的公开不可达性、Vector Greedy factory 错配、TODO/空实现，
  没有被宣称为公开能力；
- Qt、VTK、GUI、`geo_builder`、debug visualization 始终排除。

## T32 验证结果

下列命令均在当前工作树实测，而不是沿用旧任务结论：

- reference HEAD/clean status 与 28 文件 SHA-256/631 条符号审计：通过；
- 冻结 probe 内部双跑、fixture 外部双生成和 `cmp`：26/26 任务通过；
- `cargo fmt --all -- --check`：通过；
- `cargo clippy --all-targets --all-features -- -D warnings`：通过；
- `cargo test --all-targets --all-features`：209 passed、0 failed、0 ignored；
- `cargo test --release --test parity` 与 release smoke：通过；
- `cargo test --doc --all-features`：2 passed；
- `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --all-features`：通过；
- Draft 2020-12 schema 校验：正式 fixture 与正例通过，两个反例被拒绝；
- `python3 tools/audit_pure_rust.py` 与 6 个护栏单元测试：通过；
- `cargo metadata --locked --all-features` / `cargo tree --all-features`：仅零依赖 `georbf`，
  无 `custom-build`、`links` 或 native dependency；
- `cargo package --list --locked --allow-dirty` 和 `cargo package --allow-dirty --locked`：通过，
  fixture/Rust tests 被包含，仓库报告保留为开发文档，C++/reference/oracle 未进入包；
- `git diff --check`：通过。

因此 T32 可以完成并进入唯一后继 T33；性能仍为未判断。
