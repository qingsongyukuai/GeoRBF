# T34 最终发布验收

## 最终结论

T34 全部发布门槛通过。验收对象基于冻结 Surfe 提交
`290dbe0ab344f4258a4935f05cad0f153f0f69a4`；reference HEAD 精确匹配且工作树干净。
GeoRBF 已具备冻结非可视化核心的行为 parity、零生产依赖纯 Rust 边界、不低于 Surfe
的固定 workload 性能证据，以及 Linux、macOS、Windows 三平台发布 CI 证据。

## 四项最终门槛

| 状态字段 | 可复核证据 | 状态 |
|---|---|---|
| `behavior_parity` | [parity-report.md](parity-report.md)；26 个冻结 probe 的正式 fixture、74 个 family、214 个完整测试及 15 个 release parity/performance 测试全部通过。 | PASS |
| `pure_rust` | 零依赖 Cargo graph、native/build-script guard、95-file 发布包审计与从包内试构建；三平台 CI 重复执行。 | PASS |
| `performance_not_lower_than_surfe` | [performance-report.md](performance-report.md)；T34 三轮复核 12 组全部 PASS，最差 GeoRBF/Surfe ratio `0.910`。 | PASS |
| `release_audit` | 本报告、[source-traceability.md](source-traceability.md)、许可证/NOTICE、用户文档、发布包和三平台 CI。 | PASS |

## 发布物与许可证检查清单

- 根 `LICENSE`：GeoRBF 自身 2026 MIT。
- 根 `NOTICE`：包含冻结 Surfe 2017 Government of Canada MIT notice、固定提交与排除项。
- 根 `README.md`：记录安全 Builder/FittedModel 生命周期、五模型、错误边界、兼容说明和复现入口。
- [source-traceability.md](source-traceability.md)：每个生产 Rust 模块恰有一项精确
  `path@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 来源映射。
- [licensing-and-rust-boundary.md](licensing-and-rust-boundary.md)：记录 Surfe MIT 许可证继承、
  CMake/Eigen/Qt/VTK/pybind 边界与纯 Rust 发布约束。
- Cargo package 包含法律、用户、parity、性能与审计文档；不含 C++、reference、oracle、
  Eigen、Qt、VTK、CMake 产物或 native library。

## 平台矩阵

GitHub Actions [CI run 31717949048](https://github.com/qingsongyukuai/GeoRBF/actions/runs/31717949048)
在提交 `eee2c841c055f0dbaaca7fd6d77168d5fb1bf34b` 上完成，结论为 `success`。
三个 job 均执行并通过 fmt、严格 Clippy、全部 tests/doc tests、release build、严格 rustdoc、
release/Python guards、纯 Rust 审计、metadata/package 审计、dependency tree 与 package build。

| 平台 | CI 证据 | 状态 |
|---|---|---|
| Linux | [job 94507340105](https://github.com/qingsongyukuai/GeoRBF/actions/runs/31717949048/job/94507340105) | PASS |
| macOS | [job 94507340508](https://github.com/qingsongyukuai/GeoRBF/actions/runs/31717949048/job/94507340508) | PASS |
| Windows | [job 94507339687](https://github.com/qingsongyukuai/GeoRBF/actions/runs/31717949048/job/94507339687) | PASS |

`.github/workflows/ci.yml` 的 SHA-256 为
`72ae753999491d10e8388d34ba0423740f4a24eb3524788c29dac9e2d90d74fe`，矩阵精确包含
`ubuntu-latest`、`macos-latest`、`windows-latest`。

## T34 验证记录

本地实测环境为 Linux x86_64、`rustc/cargo 1.97.1`、Python 3.12.3：

- `cargo fmt --all -- --check`、严格 Clippy、release build 和 `-D warnings` rustdoc 通过。
- `cargo test --all-targets --all-features`：214 passed、0 failed、0 ignored；
  `cargo test --doc --all-features`：2 passed；release parity/performance harness：15 passed。
- Python pure-Rust/performance/release guards：16 passed；`tools/audit_pure_rust.py` 输出
  `pure-Rust audit passed`。
- Cargo metadata/tree 仅 1 个 package/node、0 dependencies、`links=null`、无
  `custom-build`；tracked native scan 为 0。
- `cargo package --allow-dirty --locked` 打包并从包内构建成功：95 files；报告自身位于包内，
  因此不记录会因证据更新而自引用漂移的 archive size/hash。
- Draft 2020-12 schema：正式 fixture/正例通过，两个反例被拒绝。26 个冻结 probe
  `--verify-repeat` 再生成与正式 fixture 逐字节相同，SHA-256 均为
  `6c26cb6bdd1f2bdba00f06c9dcdb27415cc45ecb1fb4d4e1569cc1fc3a754c28`。
- T33 comparator 三轮交替聚合全部 12 组 PASS；1-thread 六阶段 ratio 为
  `0.465/0.667/0.910/0.336/0.373/0.503`，2-thread evaluation/end-to-end ratio 为
  `0.436/0.380/0.510`。
- `cargo tree --all-features`、`cargo package --list --locked --allow-dirty`、
  `cargo package --allow-dirty --locked`、native/build-script 二次审计和 `git diff --check`
  全部通过。

三平台真实 runner 结果与本地复验共同关闭平台差异、发布包和最终四项 definition of done；
没有缺失证据、未运行门槛或未解释的 release-blocking mismatch。
