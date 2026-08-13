# T34 最终发布验收

## 当前结论

T34 的本地发布检查、T32 parity 重跑、T33 benchmark 复核和两个非 host target 的
all-target type-check 已通过；冻结来源身份为
`290dbe0ab344f4258a4935f05cad0f153f0f69a4`。实际 Linux/macOS/Windows CI 尚未运行：
远端 `origin/main` 仍为初始提交 `53766e843dbb8adc0604a11450b5ea22eaca18d5`，远端没有
当前 workflow，`gh run list --workflow CI` 报告找不到 workflow。2026-08-13 本次恢复时
`gh auth status` 已显示当前账号具备 `repo` 和 `workflow` scope，但 surfe-port 控制器明确
禁止 push，因此本会话不能发布 workflow 或触发远端矩阵。远端没有变化；T34 仍为
PENDING，`STATE.json` 不得更新为 COMPLETE。

## 四项最终门槛

| 状态字段 | 当前证据 | 状态 |
|---|---|---|
| `behavior_parity` | [parity-report.md](parity-report.md)；26 个冻结 probe 双跑再生成 fixture，与正式 fixture 字节一致；release parity 10/10。 | 本地 PASS，待三平台 CI |
| `pure_rust` | 零依赖 Cargo graph、native/build-script guard、95-file 发布包内容与试构建。 | 本地 PASS，待三平台 CI |
| `performance_not_lower_than_surfe` | [performance-report.md](performance-report.md)；T34 三轮复核全部 PASS，最差 GeoRBF/Surfe ratio `0.910`。 | PASS |
| `release_audit` | 本报告、[source-traceability.md](source-traceability.md)、许可证/NOTICE、用户文档和平台矩阵。 | PENDING：实际 CI 缺失 |

## 发布物与许可证检查清单

- 根 `LICENSE`：GeoRBF 自身 2026 MIT。
- 根 `NOTICE`：冻结 Surfe 2017 Government of Canada MIT notice、固定提交与排除项。
- 根 `README.md`：安全 Builder/FittedModel 生命周期、五模型、错误、兼容和复现入口。
- `source-traceability.md`：当前每个生产 Rust 模块恰有一项精确 `path@commit` 映射。
- Cargo package：必须包含上述法律/用户/证据文档，不得包含 C++、reference、oracle、
  Eigen、Qt、VTK、CMake 产物或 native library。

## 平台矩阵

`.github/workflows/ci.yml` 的 Linux、macOS、Windows job 对三平台运行相同的 fmt、严格
Clippy、全部 tests/doc tests、release build、严格 rustdoc、纯 Rust guard、release
audit、dependency tree 和 package build。实际 CI 运行证据将在验收后记录于此。

| 平台 | CI 证据 | 状态 |
|---|---|---|
| Linux | 本机完整 suite 通过；远端 GitHub Actions 未运行。 | PENDING CI |
| macOS | `cargo check --all-targets --all-features --target x86_64-apple-darwin` 通过；未在 macOS runner 执行测试。 | PENDING CI |
| Windows | `cargo check --all-targets --all-features --target x86_64-pc-windows-gnu` 通过；未在 Windows runner 执行测试。 | PENDING CI |

## T34 验证记录

本地实测环境为 Linux x86_64、`rustc/cargo 1.85.0`、Python 3.12.3：

- reference HEAD 精确为固定提交，工作树干净；Eigen gitlink/HEAD 为
  `36b95962756c1fce8e29b1f8bc45967f30773c00`。
- `cargo fmt --all -- --check`、严格 Clippy、release build 和 `-D warnings` rustdoc 通过。
- `cargo test --all-targets --all-features`：214 passed、0 failed、0 ignored；
  `cargo test --doc`：2 passed；release parity/performance harness：15 passed。
- Python pure-Rust/performance/release guards：16 passed；实际 `tools/audit_pure_rust.py`
  输出 `pure-Rust audit passed`。
- Cargo metadata/tree 仅 1 个 package/node、0 dependencies、`links=null`、无
  `custom-build`；tracked native scan 为 0。
- `cargo package --allow-dirty --locked` 打包并从包构建成功：95 files；报告自身位于包内，
  因此不记录会自引用漂移的 archive size/hash。
- Draft 2020-12 schema：正式 fixture/正例通过，两个反例被拒绝。26 个冻结 probe
  `--verify-repeat` 再生成与正式 fixture `cmp` 相同，SHA-256 均为
  `6c26cb6bdd1f2bdba00f06c9dcdb27415cc45ecb1fb4d4e1569cc1fc3a754c28`。
- T33 comparator 三轮交替聚合全部 12 组 PASS；1-thread 六阶段 ratio 为
  `0.465/0.667/0.910/0.336/0.373/0.503`，2-thread evaluation/end-to-end ratio 为
  `0.436/0.380/0.510`。
- CI YAML 可解析且矩阵精确包含 `ubuntu-latest`、`macos-latest`、
  `windows-latest`；workflow SHA-256 为
  `72ae753999491d10e8388d34ba0423740f4a24eb3524788c29dac9e2d90d74fe`。这只证明配置，
  不冒充未发生的 runner 结果。

因平台 CI 证据缺失，当前没有更新四项 final definition of done，也没有创建 T34 完成提交。
恢复条件是在本控制器之外把当前代码与 workflow 发布到远端并取得 Linux/macOS/Windows
矩阵结果；之后重新调用 `$surfe-port continue` 只读复核运行证据。
