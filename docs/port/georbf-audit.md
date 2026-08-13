# T05 GeoRBF 现状与纯 Rust 工程护栏审计

## 身份、范围与结论

本审计只处理 T05，Surfe 来源固定为
`290dbe0ab344f4258a4935f05cad0f153f0f69a4`。reference 按既定优先级解析到
被忽略的 `.cache/surfe-reference`；审计时 HEAD 精确匹配固定提交且工作树干净。

GeoRBF 在 T05 开始前没有 `Cargo.toml`、`Cargo.lock`、`src/` 或 CI。仓库已有的
Rust 资产只有 T03/T04 的零依赖 oracle envelope 与 golden-fixture parity helper，
因此没有需要替换的算法、公共 API 或 Cargo 依赖。T05 复用这两个 helper，建立一个
零生产依赖、零 dev/build 依赖的 library crate，并只增加工程护栏；没有实现 T06 以后
的参数、约束、核、矩阵、求解器或模型行为。

## 冻结 Surfe 构建与类型边界证据

本任务实际阅读并核对了以下冻结来源：

- `.gitmodules@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：gitlink 只有
  `pybind11` 和 `eigen-git-mirror`。
- `CMakeLists.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`：设置 C++11，
  `find_package(OpenMP)`、`find_package(Eigen3 REQUIRED)`，构建共享
  `math_lib`/`surfe_lib` 和 `surfepy`；只有 `GEO_BUILDER` 分支查找 VTK、Qt5 并
  构建 `geo_builder`/可视化 test。冻结提交不存在 `surfe_lib/CMakeLists.txt`，因此
  没有把计划中的条件性路径伪记为已读文件。
- `surfe_lib/surfe_lib_module.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
  的 `SURFE_LIB_EXPORT`、`SURFE_LIB_NO_EXPORT`、`SURFE_LIB_DEPRECATED*`，以及
  `surfe_lib/surfe_api.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的
  `class SURFE_LIB_EXPORT Surfe_API`：这是 Windows DLL/C++ ABI 边界，不迁移成
  Rust FFI。
- `math_lib/math_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
  的 `MatrixXd`/`VectorXd` QP 状态、KKT 矩阵和 `partialPivLu().solve`；
  `surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
  的 `System_Solver`、`Linear_LU_decomposition`、
  `Quadratic_Predictor_Corrector`、`Quadratic_Predictor_Corrector_LOQO`。
- `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
  的 `Matrix3f _Transform`、`SelfAdjointEigenSolver<Matrix3f>`、
  `Lagrangian_Polynomial_Basis` 的 `VectorXd`/`MatrixXd`；
  `surfe_lib/modeling_methods.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
  的矩阵抽象入口。
- `surfe_lib/{single_surface,lajaunie,stratigraphic_surfaces,continuous_property,vector_field}.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
  的 matrix/vector 参数和局部装配对象；
  `surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
  的 `MatrixXd`/`VectorXd`/`Vector3d` 公共适配；
  `surfe_pybindings/pybindings.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
  的 Eigen/pybind11 适配入口。

这些使用点按职责归入后续固定任务；T05 只把 Eigen、CMake、OpenMP、pybind11、
Qt、VTK 和 C++ ABI 判定为“由纯 Rust 实现替代”或“可视化排除”，不迁移构建系统，
也不选择后续数值算法。

## 现有 GeoRBF 逐项复用判定

| T05 前文件/资产 | 判定 | T05 处理与证据 |
| --- | --- | --- |
| `LICENSE` | 原样复用 | 项目 MIT 许可证进入 crate 包；上游 notice 的最终布局仍由 T34 审计。 |
| `tests/common/oracle_protocol.rs` | 原样复用并接入 | 保持标准库实现和“不发现/启动 oracle”边界，由 `tests/protocols.rs` 纳入 Cargo。 |
| `tests/common/parity/mod.rs` | 复用并最小扩展 | 只更新已存在 Cargo 工程的模块说明，由 `tests/protocols.rs` 纳入 Cargo；比较语义未改。 |
| `tests/fixtures/schema/*` | 原样复用 | 作为合成 schema 正反例进入测试和 crate 包；仍不是正式数值 golden。 |
| `docs/port/*` | 原样保留 | 继续作为迁移事实、协议和状态来源，不编译进生产库。 |
| `.cache/surfe-reference`、`.cache/surfe-oracle` | 保持 ignored | 不进入 Cargo metadata、正常测试或发布包。 |
| Cargo/`src`/CI | 原先不存在 | 新建最小工程边界；无可兼容实现可替换，也未推倒任何已有 API。 |

## 固定目标模块图

下图只记录 `PLAN.md` 已固定的依赖方向，不表示 T05 已提前实现任何模块：

```text
安全公共 API（T29–T30）
  -> 五模型与实际可达 Greedy（T22–T28、T31）
     -> layout / functional / assembly（T15–T17）
        -> LU / QP / LOQO / reconstruction（T18–T21）
           -> kernel / polynomial（T10–T14）
              -> parameters / errors / constraints / spatial（T06–T09）

全局 fixture parity（T32） -> 保持 parity 的性能（T33） -> 发布验收（T34）
```

当前 `src/lib.rs` 只声明 crate 文档和 `unsafe_code = "forbid"` 边界。未来模块必须按
上述顺序增加并通过同一护栏，不能用 FFI 或 native build 途径绕过纯 Rust 实现。

## Cargo 与最小依赖选择

- 单 package library crate：`georbf`，edition 2021，MSRV `1.82`。现有 parity
  helper 使用 `Option::is_none_or`/`Result::is_ok_and`，因此 `1.82` 是不重写其语义
  的明确下界；本任务在 `rustc/cargo 1.85.0` 实测。
- `[dependencies]` 为空，也没有 dev/build dependency、feature、`build.rs` 或 Cargo
  `links`。T05 的协议与 JSON helper 已可只用标准库完成，没有证据支持引入 crate。
- 后续若确需纯 Rust crate，必须记录职责与许可证，并让 all-feature metadata、全部
  target、依赖源码和发布包重新通过护栏；默认关闭 feature 不能隐藏不合规依赖。
- `Cargo.lock` 纳入版本控制，CI 与审计使用 `--locked`，使依赖图变化成为显式 diff。

## 本地和 CI 机器护栏

`tools/audit_pure_rust.py` 只使用 Python 标准库，执行以下互补检查：

1. 解析 `cargo metadata --locked --all-features --format-version 1` 的完整 resolve graph，
   拒绝已知 C/C++ compiler helper、CMake、bindgen、cxx、Eigen、Qt、VTK、pybind11、
   BLAS、OpenBLAS、LAPACK 和 MKL package family。
2. 对每个直接/传递 package 拒绝 `custom-build` target 和任何 Cargo `links`，并扫描
   package Rust 源码中的 native `extern` ABI 与 `#[link]`；未知名称的 native wrapper
   不能只靠改名绕过。
3. 检查全部 Cargo manifest 的显式 `links`，检查 tracked/untracked 非 ignored 路径，
   拒绝 `build.rs`、C/C++ 源头、CMake 输入、对象/库/可执行制品和 reference/oracle。
4. 解析 `cargo package --list --locked --allow-dirty`，对实际发布包重复路径审计。

`tests/test_pure_rust_guard.py` 用合成 metadata、源码和路径验证所有冻结禁止类别、
`custom-build`、`links`、native FFI、native 制品及 reference/oracle 均会失败，同时用
`success`、`accounting`、`black-box`、`quote`、`vtkio` 防止无边界的子串误报，并对
当前仓库运行完整护栏。

`.github/workflows/ci.yml` 在 `ubuntu-latest`、`macos-latest`、`windows-latest` 上固定
运行 stable Rust 的 fmt、严格 clippy、全部 target/feature tests、doc tests、Python
护栏测试、实际护栏、依赖树和发布包清单。workflow 不 checkout reference，也没有
oracle discovery、CMake、C/C++ 或 native 工具链步骤。T05 只交付并本地验证 workflow
所执行的命令，不把尚未发生的远端平台运行伪记为 CI 通过；T34 仍须审计实际平台证据。

## 发布包内容基线

T05 的 `Cargo.toml` 使用显式 `include` allowlist。`cargo package --list --locked
--allow-dirty` 的基线为：

```text
.cargo_vcs_info.json
Cargo.lock
Cargo.toml
Cargo.toml.orig
LICENSE
src/lib.rs
tests/common/oracle_protocol.rs
tests/common/parity/mod.rs
tests/fixtures/schema/golden-fixture-v1.schema.json
tests/fixtures/schema/invalid-ambiguous-tolerance.json
tests/fixtures/schema/invalid-wrong-source.json
tests/fixtures/schema/valid-minimal.json
tests/protocols.rs
```

清单不含 C/C++、Eigen、Qt、VTK、CMake、build script、reference、oracle、Python 审计
工具、CI 配置或迁移控制文档。审计工具与 CI 是仓库级发布门槛，而不是生产 crate
运行时内容。提交后还必须在干净工作树以不带 `--allow-dirty` 的正式命令复核此清单。

## T05 验证命令

```text
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
cargo test --doc --all-features
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -v tests/test_pure_rust_guard.py
python3 tools/audit_pure_rust.py
cargo tree --all-features
cargo package --list --locked --allow-dirty
git diff --check
```

提交前结果必须全部通过；干净提交后再执行计划要求的精确
`cargo package --list --locked` 并记录最终结果。T05 不执行数值 parity 或性能比较，
也不把这两项写成通过。
