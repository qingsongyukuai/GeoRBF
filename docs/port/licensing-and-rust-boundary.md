# 许可证、来源追踪与纯 Rust 边界

## 上游许可证结论

冻结的 `License.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 声明：除另有注明
外，SURFE 程序源码属于加拿大政府 Crown Copyright，并按 MIT License 分发；其
版权声明为 `Copyright (c) 2017 Government of Canada`。MIT 条款允许使用、复制、
修改、合并、发布、分发、再许可和销售，但任何源码副本或实质性部分必须保留该版权
声明和许可声明，且软件不附带担保。

上游还明确指出 Canada wordmark/相关图形不获授权复用，第三方软件权利归各自权利人。
GeoRBF 不迁移 wordmark、相关图形、上游截图或第三方 submodule 内容。

GeoRBF 根 `LICENSE` 当前是项目自身的 MIT 文本和 2026 项目版权声明。T00 不把上游
源码或实质性源码副本加入仓库，因此不改写根许可证。只要后续提交包含翻译、改编或
其他可能构成上游源码实质性部分的内容，就必须同时保留上游 2017 Government of
Canada MIT notice；最终 NOTICE/许可证布局和发布包复核由 T34 完成。不得以“纯 Rust
重写”为理由删除适用的上游 notice。

本结论是仓库内迁移规则，不替代针对发布场景的专业法律意见。

## 来源追踪规则

每个迁移任务必须在 `docs/port/JOURNAL.md` 记录：

- 精确的 `path@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；
- 实际阅读和迁移的类型/函数/分支；
- 对应 GeoRBF 文件；
- oracle/fixture 生成身份与验证命令；
- 与上游不同的安全处理及其兼容证据。

新增 Rust 模块若直接翻译或紧密改编冻结实现，应在模块文档或邻近来源清单使用同样
的 path、commit 和 symbol 记录。不能只记录仓库 URL，也不能把数学论文或 README
误写成实际源码来源。T34 汇总逐模块来源清单并检查上游 notice 是否随源码和 crate
发布。

正式 fixture 只能包含重现行为所需的数值/离散数据、生成元数据和提交身份，不得
嵌入上游 C++ 源码。fixture 许可和来源在 T04/T32 单独审阅。

## 第三方边界

冻结 `.gitmodules@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 声明 Eigen mirror
和 pybind11；冻结 `CMakeLists.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
使用 Eigen3、pybind11、可选 OpenMP，并在 `GEO_BUILDER` 分支使用 VTK 与 Qt5。
这些依赖的存在只说明 Surfe 构建关系，不授权复制，也不允许它们进入 GeoRBF 生产
依赖。GeoRBF 不 vendor 或迁移这些 submodule/库。

任何新增纯 Rust crate 仍须单独审阅许可证、`build.rs`、Cargo `links` 元数据和完整
依赖树；“crate 是 Rust API”不等于没有 native 代码。T05 建立机器可检护栏，T34
在发布包上复核。

## 生产、测试与发布禁止项

GeoRBF 的正常 debug/release 构建、全部 Cargo 测试、文档测试和发布包不得依赖或
编译以下任何项：

- C 或 C++ 生产源码及 C++ ABI；
- Eigen；
- Qt；
- VTK；
- CMake 构建或运行时依赖；
- bindgen；
- `cxx`/CXX bridge；
- 用于生产的 `cc` build dependency；
- MKL；
- OpenBLAS；
- LAPACK；
- native BLAS FFI；
- pybind11、冻结 Surfe reference 或 C++ oracle。

禁止项适用于直接、可选、feature-gated、build、dev 和传递依赖：不能靠默认关闭
feature 隐藏。正常 `cargo test --all-targets --all-features` 和打包试构建必须仍为
纯 Rust。纯 Rust 实现可以使用标准库或经审计、没有 native/build-script 边界的
Rust 依赖。

reference/oracle 唯一允许的位置为外部 `SURFE_REFERENCE_DIR`、`../surfe`，或被
忽略的 `.cache/surfe-reference`/`.cache/surfe-oracle`。它们仅允许被显式 oracle、
fixture 生成和 benchmark 命令调用，不能被 `build.rs`、Cargo feature、测试自动发现
逻辑或发布脚本设为成功前置条件。

## 发布包与审计判定

T05 起至少检查：

- `cargo tree --all-features` 中的直接、传递和 build dependencies；
- 所有 `Cargo.toml`、`Cargo.lock`、`build.rs`、Cargo `links` 和 native 编译命令；
- `cargo package --list` 与实际打包试构建；
- 仓库中被跟踪的 C/C++、对象、静态/共享库、reference/oracle 路径；
- Linux、macOS、Windows 支持矩阵上的构建和测试。

只有机器检查确认禁止项不存在、正常测试不需要 reference、发布包不含 C++/oracle，
才可把纯 Rust 门槛标为通过。T00 仓库尚无 Cargo 工程，因此只冻结边界，不声称已
完成 T05 或 T34 的依赖/发布审计。
