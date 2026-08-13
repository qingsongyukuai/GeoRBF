# Surfe Port Journal

本日志只追加已执行任务的事实与证据。任务状态以 `STATE.json` 为唯一事实来源。

## INIT — 初始化 surfe-port 持久化控制器

- 日期：2026-08-13。
- 状态：完成。
- Surfe 源码：未克隆、未构建、未建立 reference/oracle、未执行 T00；仅为固定计划的文件/符号映射核对了公开冻结提交的 `surfe_lib`、`math_lib` 和 `test` 目录及相关声明。
- 修改文件：`.agents/skills/surfe-port/SKILL.md`、`.agents/skills/surfe-port/agents/openai.yaml`、`.gitignore`、`docs/port/PLAN.md`、`docs/port/STATE.json`、`docs/port/JOURNAL.md`。
- 核心实现：无。仅建立显式触发的仓库级 Skill、T00–T34 固定计划、跨会话 JSON 状态、追加式日志和 reference/oracle 忽略规则。
- 验证命令：Skill `quick_validate.py`；JSON 解析和状态一致性检查；YAML frontmatter 与 `openai.yaml` policy 检查；T00–T34 连续性和每任务必备字段检查；`git diff --check`；变更范围检查。
- 验证结果：全部通过；`next_task` 为 `T00`，未产生算法代码、正式 golden fixture、C++ 源码或 T00 交付物。
- Parity 证据：不适用；INIT 不执行 parity。
- 性能证据：不适用；INIT 不进行性能判断。
- 后续发现：无；下一任务固定为 T00。

## T00 — 冻结 Surfe reference、许可证、迁移范围、兼容策略和纯 Rust 边界

- 日期：2026-08-13。
- 状态：完成。
- Surfe reference：从官方仓库克隆到被忽略的 `.cache/surfe-reference`；detached HEAD 与 `git rev-parse HEAD` 均确认 `290dbe0ab344f4258a4935f05cad0f153f0f69a4`，reference 工作树干净。
- 阅读的 Surfe 源码：`License.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`README.md@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`.gitmodules@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`CMakeLists.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `Surfe_API::ComputeInterpolant`、`EvaluateInterpolantAtPoint(s)`、`EvaluateVectorInterpolantAtPoint(s)` 和约束/参数/状态入口，以及 `surfe_lib/debug.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `open_console_window` 排除证据。
- 修改文件：`docs/port/source-baseline.md`、`docs/port/compatibility.md`、`docs/port/licensing-and-rust-boundary.md`、`docs/port/STATE.json`、`docs/port/JOURNAL.md`。
- 核心实现：固定上游仓库/提交与 reference 三段优先级；记录 `License.txt` SHA-256 `9fd4e80cac11aa3d00278e4d2634ff0c6a169303014155f51873f9c4e6e6441f`、Crown Copyright/MIT notice 义务和逐模块来源追踪规则；对冻结仓库顶级路径给出纳入、证据限定或排除判定；固定有效输入、数值/离散 parity、缺陷和类型化错误策略；列出生产、测试和发布的全部纯 Rust/native 禁止项。未建立 oracle、fixture 或算法实现。
- 验证命令：`git -C .cache/surfe-reference rev-parse HEAD` 与 reference clean-status 检查；工作树及 `git show <commit>:License.txt` 双重 SHA-256/许可证关键文本校验；`git check-ignore`；固定提交、reference 优先级、排除项、兼容规则和全部 native 禁止项的 `rg` 机器检查；tracked C/C++/二进制拒绝检查；T00 状态/任务序列 `jq`/`diff` 校验；`git diff --check`。
- 验证结果：全部通过；reference 提交精确匹配且无本地修改，reference/oracle 路径被忽略，三份交付文档的必需字段齐全，GeoRBF 未跟踪 Surfe 源码、对象、库或可执行文件。
- Parity 证据：T00 只完成不可漂移的 source/license/scope 判定基线；固定提交和许可证内容已逐字节验证。按计划尚未构建 T03 oracle、T04 fixture 或执行数值 parity，未将其误报为通过。
- 性能证据：不适用；T00 不进行性能判断或性能声明。
- 后续发现：`Surfe_API::SetTangentConstraints` 调用 planar 添加入口、`Surfe_API(int)` 未显式初始化全部状态布尔值、`SetRegressionSmoothing`/`SetGreedyAlgorithm` 无条件启用等源码观察已登记到既有 T02/T30/T31 核验范围；未创建任务或提前修复。下一任务固定为 T01。

## T01 — 建立完整 C++ 清单、调用链、数据流和依赖图

- 日期：2026-08-13。
- 状态：完成。
- Surfe reference：按优先级解析到被忽略的 `.cache/surfe-reference`；`git rev-parse HEAD` 再次确认 `290dbe0ab344f4258a4935f05cad0f153f0f69a4`，reference 工作树无修改。
- 阅读的 Surfe 源码：完整审计 `surfe_lib/*.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`math_lib/*.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_pybindings/pybindings.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的核心适配调用、`test/main.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`，并读取根 `CMakeLists.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 以固定 Eigen/OpenMP/pybind11 与可选 Qt/VTK/geo_builder 边界。重点逐体核对 `Surfe_API` 两个工厂和拟合/评估入口、`GRBF_Modelling_Methods::{get_method,setup_basis_functions,run_greedy_algorithm}`、五模型的全部虚函数覆写、`RBFKernel`/`Modified_Kernel` 组合、三个 solver 与 `Math_methods` QP 路径。
- 修改文件：`docs/port/inventory.md`、`docs/port/call-graph.md`、`docs/port/data-flow.md`、`docs/port/STATE.json`、`docs/port/JOURNAL.md`。
- 核心实现：无 Rust 算法。建立 30 个范围文件的逐文件清单、类型/字段/逐符号清单、五模型和 12 个纯虚槽覆写矩阵；记录 API/Greedy 双工厂、从约束到分组/pair/layout/matrix/RHS/weights/scalar/gradient 的字段级流、LU/QP/LOQO/重建路径、裸指针所有权，以及 Eigen/CMake/OpenMP/pybind11/Qt/VTK/geo_builder 的迁移/替代/排除判定。
- 验证命令：`git -C .cache/surfe-reference rev-parse HEAD` 与 clean-status；`find` 枚举范围文件并用 `rg -F` 将每个 `path@commit` 与文档双向核对；`rg` 提取 class/struct/enum；Perl 去注释后提取限定符号并以 `rg` 核对清单；表驱动检查五模型在公开工厂中的分支和每个模型头中的 12 个覆写；检查三份交付物包含调用者、被调用者、状态所有者、全部 layout、依赖处置和后续任务归属；`git diff --check` 与 staged 变更范围检查。
- 验证结果：全部 T01 门槛通过；28 个核心 `.h/.cpp` 加 pybindings/test 两个范围文件均有准确提交路径，源码类型、枚举、去注释限定符号、五个公开工厂分支和 60 个模型覆写槽核对无缺口；三份文档明确入口、调用关系、字段所有者及每项依赖的迁移/Rust 替代/排除判定。额外的非门槛 C++ 编译探针因 reference 的 Eigen submodule 未初始化而未运行，未把它记录为构建通过；T01 的固定验收不要求构建，T03 oracle 构建仍必须独立取得实测证据。
- Parity 证据：T01 的离散 parity 证据仅为冻结声明/定义/调用点与清单的集合核对；未建立 T03 oracle、T04 fixture 或数值 golden，未宣称数值 parity。
- 性能证据：不适用；T01 不实现或基准测试算法，未作性能判断。
- 后续发现：公开 API 工厂为五个显式分支，而 `GRBF_Modelling_Methods::get_method` 的 Greedy 工厂无显式 Vector Field 分支；若干内联常量体、抛出体、仅声明符号、`SetTangentConstraints` 的实际调用目标、批量状态检查/共享可变核与状态位初始化等事实已映射到既有 T02/T30/T31，未提前分类或修复。下一任务固定为 T02。

## T02 — 逐符号能力分类

- 日期：2026-08-13。
- 状态：完成。
- Surfe reference：按优先级解析到被忽略的 `.cache/surfe-reference`；`git rev-parse HEAD` 精确为 `290dbe0ab344f4258a4935f05cad0f153f0f69a4`，reference 工作树干净。
- 阅读的 Surfe 源码：完整复核 `surfe_lib/*.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`math_lib/*.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_pybindings/pybindings.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 和 `test/main.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的声明、定义与活动调用点。重点逐体核对 `basis.cpp` 的 15 类核方法、`modeling_methods.cpp` 的 `GRBF_Modelling_Methods::{check_input_data,run_greedy_algorithm}`、五模型 60 个覆写槽、三个完整 `convert_modified_kernel_to_rbf_kernel`、`continuous_property.cpp` 全部覆写、`surfe_api.cpp` 的状态/setter/批量入口、`modelling_input.cpp` 的空输入与 Greedy residual selector，以及 `grbf_exceptions.h` 的 23 个具体错误类别和 `SurfeExceptions` wrapper。
- 修改文件：`docs/port/symbol-classification.md`、`docs/port/compatibility.md`、`docs/port/STATE.json`、`docs/port/JOURNAL.md`。
- 核心实现：无 Rust 或 C++ 算法实现。为 T01 清单建立互斥的 `I/P/T/U/D/V/X` 逐符号分类，每项记录定义/调用证据、可观察缺口和唯一既有任务归属；把源级缺陷同步到兼容策略，未修改算法、未新增任务。
- 验证命令：`git -C .cache/surfe-reference rev-parse HEAD/status --short`；`find` 核对 30 个范围文件；Perl 去注释限定符号提取与 `rg` 反向检索；从 T01 清单提取代码标识符并反向检索分类表；表驱动核对 13 个完整核的 208 个 `.cpp` 公式和 13 个内联 clone、`R/AR` 的 30 个哨兵导数符号、五模型 60 个覆写槽；`rg` 搜索 TODO/空体/恒定返回/throw/调用点；`jq` 状态序列校验；`git diff --check`。
- 验证结果：全部通过。源级集合核对覆盖 501 个 T01 代码标识符 token 和 432 个去注释限定定义/调用 token；完整核共 221 个公式/clone 符号，`R/AR` 为 24 个字面 `throw -666` 体加 6 个间接别名；`run_greedy_algorithm` 只有声明/定义文件，三个完整 reconstruction 函数零调用，`use_greedy` 的五处命中均是声明、入参、默认值或写入，没有消费点。
- Parity 证据：T02 只固定源码能力/不可达/缺陷分类；未建立或运行 T03 C++ oracle，未生成 T04 golden，未宣称数值 parity。
- 性能证据：不适用；T02 不实现或基准测试算法，未作性能判断。
- 后续发现：`Parameters::use_greedy` 和 `GRBF_Modelling_Methods::run_greedy_algorithm` 不可达归 T31；Single/Lajaunie/Stratigraphic 的完整 reconstruction 函数不可达归 T21；状态、约束路由、越界和共享可变 kernel 问题分别归既有 T26/T27/T29/T30/T31。T01 文本的“22 个异常”是计数笔误，实际列出 23 个；计划提示的 `Point::operator==` 在冻结声明/定义中不存在，真实同点符号是 `collocated`。以上均映射到既有任务或清单更正，未创建新任务。下一任务固定为 T03。

## T03 — 建立仓库外冻结 Surfe reference oracle

- 日期：2026-08-13。
- 状态：完成。
- Surfe reference：按既定优先级解析到被忽略的 `.cache/surfe-reference`；`git rev-parse HEAD` 精确为 `290dbe0ab344f4258a4935f05cad0f153f0f69a4`，reference 工作树干净。为外部构建初始化冻结 gitlink `eigen-git-mirror@36b95962756c1fce8e29b1f8bc45967f30773c00`；首次网络传输中断后 Git 自动重试成功，最终 submodule HEAD 与 gitlink 精确一致。
- 阅读的 Surfe 源码：`surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的约束 setter、`ComputeInterpolant`、scalar/vector 单点与批量入口及错误状态；`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `RBFKernel::{basis,dx_p1,dx_p2,dy_p1,dy_p2,dz_p1,dz_p2,dxx..dzz,get_global_anisotropy,scaled_radius,basis_*}`、`Lagrangian_Polynomial_Basis` 和 `Modified_Kernel::basis_*`；`surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 LU/QP/LOQO 输入、验证、`solve` 和 `weights`；`surfe_lib/modeling_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的清洗、level/group/reference、basis、equality matrix 和所有权；`surfe_lib/{single_surface,lajaunie,stratigraphic_surfaces,continuous_property,vector_field}.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `process_input_data/get_method_parameters/get_interpolation_matrix/get_*values/setup_system_solver/eval_*` 及 increment/层序提取点；必要依赖 `surfe_lib/modelling_input.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`math_lib/math_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；以及 `test/main.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的调用顺序和可视化排除证据。
- 修改文件：`docs/port/oracle-protocol.md`、`docs/port/oracle-manifest.json`、`tests/common/oracle_protocol.rs`、`docs/port/STATE.json`、`docs/port/JOURNAL.md`。外部 `.cache/surfe-oracle/oracle_smoke.cpp`、request/response、Rust test binary 和 C++ executable 均被忽略且未提交。
- 核心实现：固定 `georbf-surfe-oracle` v1 UTF-8 JSON Lines envelope、commit 身份、成功/错误不变量、binary64/binary32/non-finite 编码、row-major matrix/vector 表示、稳定错误结构和确定性环境；定义 `identity`、`kernel.evaluate`、`model.run`、`solver.run`、`error.probe`，覆盖核值、两点一阶导数、混合 Hessian、anisotropy、Modified Kernel、五模型约束/分组/layout/matrix/RHS、LU/QP/LOQO 求解证据、scalar/gradient 和错误。新增不发现或启动 oracle 的纯标准库 Rust envelope 校验辅助；正常构建仍不依赖 reference。
- 外部构建与 smoke：使用 `g++ 13.3.0`、`-std=c++11 -O2 -DNDEBUG -fopenmp` 和冻结 Eigen，编译冻结 `math_lib/surfe_lib` 与被忽略 adapter。`model.run/single_surface_linear_smoke` 经 `Surfe_API::ComputeInterpolant` 运行四个 interface、一个 planar、Cubic、一阶 polynomial，并与直接模型矩阵证据路径交叉核对；输出 11×11 matrix、11 项 RHS/weights、LU 残差和 scalar/gradient。adapter/source、请求、响应和 executable 的最终 SHA-256 分别为 `4e3f56fabd073f7ed489cceaf3f04ee4df30dbb9169207ec7767efa76ffd956a`、`cd1a994c0d075f3cc001c5a1518b4b691eb0eca13bcbbc0370385b26043cb267`、`b15aaba2405d51278d61d1d4b0d5e57c0c4fcb2ae1adc1cecabfbd067321cc7c`、`cd9815aaf26a60d725883bb2a1dff88f24c313ef1b94ce1a9685cece6c399d11`。
- 验证命令：reference HEAD/clean-status 与 Eigen gitlink/HEAD 检查；上述 `g++` 构建；固定 `OMP_NUM_THREADS=1 LC_ALL=C TZ=UTC` 连续运行相同请求并 `cmp`；`jq` 校验 request/response envelope、commit、成功/错误 schema、11×11/121 项 matrix、11 项 RHS/weights、有限性、solver 分支与 prediction；`rustc --test tests/common/oracle_protocol.rs` 并把真实 smoke 行显式传入；manifest `jq` 覆盖检查；`git check-ignore`、tracked native 扩展拒绝检查、`git status --ignored`；状态连续性检查；`git diff --check`。
- 验证结果：全部 T03 门槛通过。两次成功响应逐字节相同且 SHA-256 均为 `b15aaba2405d51278d61d1d4b0d5e57c0c4fcb2ae1adc1cecabfbd067321cc7c`；response commit 精确匹配，所有数值有限，LU `attempted/success` 为真，`residual_l2 = 0`，公开 API prediction 为 scalar `0.5`、gradient `[0,0,1]`；实际 invalid request 返回 request-stage `invalid_request`；Rust 4 个协议测试全部通过；仓库未跟踪任何 C/C++/Eigen/object/library/executable/response。
- Parity 证据：T03 证明冻结 oracle 可构建、公开 API 与中间 matrix/RHS/solve/prediction 可在同一协议响应中交叉提取，且同输入输出字节稳定；它不生成正式 golden、不定义 T04 容差，也不宣称 T32 全局 parity。smoke 的公开 API 与直接模型路径在 `1e-14` 内一致，实测结果恰好相同。
- 性能证据：不适用；T03 只验证 oracle 可运行和确定性，未执行同机 Surfe/GeoRBF benchmark，未作性能声明。
- 后续发现：构建再次触发冻结 `modelling_input.cpp` tangent residual selector 的“non-void function control reaches end”警告，与 T02 已登记的不可达缺陷一致，仍归 T31，不新建任务、不修改 reference。下一任务固定为 T04。

## T04 — 定义 golden fixture 与差分测试协议

- 日期：2026-08-13。
- 状态：完成。恢复上下文时确认仓库根目录没有 `AGENTS.md`，按控制器规则记录后继续；任务开始前工作区干净，当前分支为 `main`。
- Surfe reference：按优先级解析到被忽略的 `.cache/surfe-reference`；`git rev-parse HEAD` 精确为 `290dbe0ab344f4258a4935f05cad0f153f0f69a4`，reference 工作树干净，Eigen gitlink 仍为 `36b95962756c1fce8e29b1f8bc45967f30773c00`。
- 阅读的 Surfe 源码：`surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `RBFKernel::{basis,dx_p1,dx_p2,dy_p1,dy_p2,dz_p1,dz_p2,dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz,get_global_anisotropy,scaled_radius,basis_*}`、九个 isotropic 核、anisotropic 核和 `Modified_Kernel::basis_*`；`surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `System_Solver`、`Linear_LU_decomposition::solve`、`Quadratic_Predictor_Corrector::solve`、`Quadratic_Predictor_Corrector_LOQO::solve` 及验证入口；`surfe_lib/{single_surface,lajaunie,stratigraphic_surfaces,continuous_property,vector_field}.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `get_interpolation_matrix`、`get_equality_values`、可达 `get_inequality_matrix/get_inequality_values`、`setup_system_solver`、`eval_scalar_interpolant_at_point` 和 `eval_vector_interpolant_at_point`。
- 修改文件：`docs/port/fixtures.md`、`tests/fixtures/schema/golden-fixture-v1.schema.json`、三个 `tests/fixtures/schema/*` 合成正反例、`tests/common/parity/mod.rs`、`docs/port/STATE.json`、`docs/port/JOURNAL.md`。
- 核心实现：定义 `georbf-surfe-golden` v1 envelope、冻结 source/oracle 身份、生成 argv/环境/平台/精度和 request/response 行 SHA-256；固定默认 exact、唯一 numeric path 分类和 12 个分层 binary64/binary32 容差；固定 signed zero、tagged NaN/Inf、row-major matrix/vector 及错误编码；以精确有序坐标、方向、levels、evaluation points、显式 LU/QP/LOQO matrix/RHS 固定九核、两点一阶导数、混合 Hessian、anisotropy、Modified Kernel、五模型、matrix/RHS、scalar/gradient、objective/iteration/error 的确定性 case family；规定非唯一权重只作诊断，以机器可读 residual/feasibility ceiling 配合实际尝试、有限性、目标和 witness predictions 共同判定 pass/fail。新增零依赖 Rust JSON reader/canonical serializer、fixture validator、matrix/vector invariant、数值/路径 comparator；未实现迁移算法，未提交正式数值 golden。
- 验证命令：reference HEAD/clean-status/gitlink 检查；`rustfmt --edition 2021 --check tests/common/parity/mod.rs`；`rustc --edition=2021 -D warnings --test tests/common/parity/mod.rs` 后运行；Python `jsonschema` Draft 2020-12 schema 自校验和一个正例/两个反例；`jq -c` 与 `sha256sum` 复核 schema 正例内 request/response 行 hash；固定环境下用 T03 外部 oracle 对相同 smoke request 连续运行两次并 `cmp`；coverage/name `rg` 检查；`git check-ignore`；tracked native/正式 golden 拒绝检查；状态连续性检查；`git diff --check`。仓库尚无 `Cargo.toml`，因此计划所述条件性 `cargo test` 不适用，未伪造为运行通过。
- 验证结果：Rust 7 个测试全部通过，覆盖 schema 正反例、规范序列化往返、容差闭边界及下一 binary64 值失败、signed zero、tagged non-finite、路径唯一分类、分层阈值和 matrix/vector 长度/编码；Draft 2020-12 schema 有效，`valid-minimal.json` 通过，source 漂移和 `all_numeric` 单一宽容差反例均被拒绝；正例 hash 与内嵌值一致。冻结 oracle 两次响应逐字节相同且 SHA-256 均为 `b15aaba2405d51278d61d1d4b0d5e57c0c4fcb2ae1adc1cecabfbd067321cc7c`，source commit/status 复核通过；所有提交外 reference/oracle 文件保持 ignored。
- Parity 证据：T04 固定了后续差分的可审计格式、离散 exact 规则、分层数值判定和确定性覆盖矩阵，并以真实冻结 oracle 重跑证明同输入输出字节稳定；schema 目录中的三个 JSON 只为合成结构测试，不是正式数值 fixture。本任务未宣称任何 Rust 数学实现、模型或 T32 全局数值 parity 已通过。
- 性能证据：不适用；T04 不实现算法或运行同机 Surfe/GeoRBF benchmark，未作性能判断。
- 后续发现：T05 建立 Cargo 工程后应把现有零依赖 parity helper 纳入普通 `cargo test`，但不得让正常测试发现或启动 oracle；该工作已属于固定 T05，不新增任务。下一任务固定为 T05。

## T05 — 审查 GeoRBF 并建立纯 Rust 工程护栏

- 日期：2026-08-13。
- 状态：完成。恢复上下文时根目录没有 `AGENTS.md`；任务开始前 `main` 分支工作区干净。`STATE.json` 的 INIT–T04 连续前缀、`last_completed_task=T04`、`active_task=null`、`next_task=T05` 校验一致后才激活 T05。
- Surfe reference：按优先级解析到被忽略的 `.cache/surfe-reference`；`git rev-parse HEAD` 精确为 `290dbe0ab344f4258a4935f05cad0f153f0f69a4`，reference 工作树干净。
- 阅读的 Surfe 源码：`.gitmodules@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 Eigen/pybind11 gitlink；`CMakeLists.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 C++11、OpenMP、Eigen3、共享 `math_lib/surfe_lib`、pybind11 和条件性 Qt5/VTK/geo_builder 分支（冻结提交不存在计划中的条件性 `surfe_lib/CMakeLists.txt`）；`surfe_lib/surfe_lib_module.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 `SURFE_LIB_EXPORT/SURFE_LIB_NO_EXPORT/SURFE_LIB_DEPRECATED*` 及 `surfe_lib/surfe_api.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 DLL 导出类；`math_lib/math_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`、`surfe_lib/{basis,matrix_solver,modeling_methods,single_surface,lajaunie,stratigraphic_surfaces,continuous_property,vector_field,surfe_api}.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 和 `surfe_pybindings/pybindings.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 Eigen `Matrix3f/MatrixXd/VectorXd/Vector3d`、`SelfAdjointEigenSolver`、`partialPivLu` 及矩阵/向量 API 使用点。
- 修改文件：`Cargo.toml`、`Cargo.lock`、`src/lib.rs`、`.github/workflows/ci.yml`、`.gitignore`、`tools/audit_pure_rust.py`、`tests/protocols.rs`、`tests/test_pure_rust_guard.py`、`tests/common/parity/mod.rs`、`docs/port/georbf-audit.md`、`docs/port/STATE.json`、`docs/port/JOURNAL.md`。
- 核心实现：先为 Cargo 接入和护栏写失败测试，再建立 edition 2021、MSRV 1.82 的单 package library crate；生产/dev/build 依赖与 feature 均为空，`unsafe_code` 为 forbid，现有 oracle/parity helper 通过 `tests/protocols.rs` 进入正常 Cargo test 且仍不发现/启动 reference。新增全标准库 Python 护栏：审计 all-feature 完整 Cargo metadata graph、package family、所有 `custom-build` target、Cargo `links`、直接/传递 Rust 源码 native ABI/`#[link]`、全部 manifest、tracked/untracked 非 ignored 路径和 `cargo package --list`，拒绝 C/C++、Eigen、Qt、VTK、CMake、bindgen、cxx、pybind11、BLAS/OpenBLAS/LAPACK/MKL、build script、native 制品及 reference/oracle。CI 在 Linux/macOS/Windows 矩阵配置 fmt、严格 clippy、全部 test/doc、护栏、tree 和 package list；没有提前实现任何数学或模型功能。
- 验证命令：`rustc/cargo/rustfmt/clippy --version`；reference HEAD/clean-status 和 `rg` Eigen/export 使用点；`cargo fmt --check`；`cargo clippy --all-targets --all-features -- -D warnings`；`cargo test --all-targets --all-features`；`cargo test --doc --all-features`；`cargo build --release --all-features`；`RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --all-features`；`PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -v tests/test_pure_rust_guard.py`；`python3 tools/audit_pure_rust.py`；`cargo metadata --locked --all-features --format-version 1`；`cargo tree --all-features`；`cargo package --list --locked --allow-dirty`；`cargo package --allow-dirty --locked`；workflow YAML parse；状态连续性与变更范围检查；`git diff --check`。提交后的干净工作树再运行精确 `cargo package --list --locked`。
- 验证结果：本地 `rustc/cargo 1.85.0` 下全部通过。严格 clippy 无 warning；Cargo 运行 11 个 protocol/fixture 测试全部通过，lib/doc test 为 0 个且成功；6 个 Python 护栏测试覆盖全部冻结禁止 family、custom build、links、native FFI、native/package/reference 路径及防误报；实际审计输出 `pure-Rust audit passed`。metadata/tree 只有 `georbf@0.1.0` 一个 package/resolve node，依赖为零、`links=null`、无 custom-build；release build 与严格 rustdoc 成功。package allowlist 为 13 个文件，不含 C/C++、reference/oracle、build script 或审计/CI 工具；试打包和从包验证构建成功，大小 68.3 KiB、压缩后 13.6 KiB。workflow YAML 可解析；没有把未发生的远端三平台 run 伪记为 CI 通过，实际平台证据仍归 T34。
- Parity 证据：本任务不实现或比较 Surfe 数学行为；T03/T04 的 4 个 oracle protocol 测试与 7 个 fixture/parity 基础测试现已在普通 `cargo test` 中全部通过，且测试不访问 oracle。未生成正式 golden，未宣称 T32 数值 parity。
- 性能证据：不适用；T05 没有算法或同机 Surfe/GeoRBF benchmark，未作性能判断。
- 后续发现：没有新增任务。T06 起新增 Rust 模块和任何经证明必要的纯 Rust 依赖都必须继续通过 T05 护栏；下一任务固定为 T06。
