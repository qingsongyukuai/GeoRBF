# GeoRBF Rust 数值后端能力边界

状态：研究结论，供后续数值政策票使用

研究票据：[调查 Rust 稠密线性、凸 QP 与 SOCP 后端的能力边界](https://github.com/qingsongyukuai/GeoRBF/issues/7)

证据快照：2026-08-01

范围：3D/f64、小中规模、最多 2,000 个 lowering 后标量约束；调查事实与风险，不选择 GeoRBF 后端、不实现产品代码

## 结论摘要

1. Rust 生态已经能覆盖本项目的三类标准问题，但不是靠一个同质后端完成：

   - `faer 0.24.4` 与 `nalgebra 0.35.0` 都有纯 Rust 的 Cholesky、Bunch–Kaufman 对称不定分解、列主元 QR 与 SVD；`faer` 更明确面向中大型高性能稠密/稀疏计算，`nalgebra` 的分解结果较易检查，但两者都没有替 GeoRBF 定义 rank/condition tolerance。
   - `nalgebra-lapack 0.28.0` 与 `ndarray-linalg 0.18.1` 能调用成熟 LAPACK；前者直接包装 `?geqp3` 并提供可选 rank heuristic，后者提供 `?sytrf` Bunch–Kaufman、SVD/`?gesdd`、`?gelsd` rank 与一般 LU 的 `?gecon` reciprocal-condition estimate。代价是 BLAS/LAPACK 原生依赖、链接策略和平台矩阵成为产品责任。
   - `Clarabel.rs 0.11.1` 在一个 Rust API 中直接覆盖 convex QP 与 SOCP，并暴露 residual、gap、infeasibility residual/certificate statuses；默认 QP/SOCP 路径可不依赖 C/Fortran。`osqp.rs 1.0.1` 是 OSQP 1.0.0 C 库的 Rust 包装，只覆盖 convex QP，但提供明确 warm start、增量数据更新和可取出的 infeasibility certificates。

2. “solver 返回成功”不足以满足 GeoRBF 契约。Clarabel 与 OSQP 都会缩放并对 KKT 求解作算法性 regularization；它们的 residual/status 使用各自标准形式和 tolerance。后续数值政策必须保留变换映射，并在 Canonical Problem IR 的物理单位中独立重算 objective、hard violations、SOC membership、CPD side conditions 与 certificates。

3. rank/condition 可见性是独立能力，不是分解名称的附赠品。列主元 QR 的对角线或 SVD singular values 能提供证据，但 threshold 仍由 GeoRBF 决定；Cholesky 失败、精确零 pivot 或 solver `NumericalError` 都不能自动改名为 `RankDeficient` 或 `Infeasible`。

4. 产品基线的 “2,000” 是 lowering 后标量约束数，不等于 $2{,}000\times2{,}000$ 后端矩阵。centers、polynomial/semantic latents、slacks、cone rows 和 homogeneous/KKT variables 都会改变实际维数；容量判断必须记录实际 primal dimension $n$、constraint/cone-row dimension $m$、KKT order $k$ 与 nonzero/fill counts。

5. 本研究不给出“选哪一个”的答案。它给下一张[数值问题标准化、后端路由与稳定性政策](https://github.com/qingsongyukuai/GeoRBF/issues/8)的直接输入：需要分别决定 dense linear baseline、rank oracle、QP/SOCP route、certificate acceptance、warm-start 是否进入 v1，以及 native-dependency/MSRV budget。

## 1. 证据规则与项目边界

### 1.1 标签

- **已验证事实**：由候选 crate 的发布文档、源码、Cargo manifest、上游 solver 文档或算法论文直接支持。
- **推论**：由已验证接口与 GeoRBF 已定契约推导；不是上游兼容或稳定性承诺。
- **未知/待测**：公开资料或本次静态检查不足以回答，必须用固定版本 build/bench/failure corpus 验证。

### 1.2 GeoRBF 已确定的输入

本研究遵守以下既有边界：

- Canonical Problem IR 是物理单位中的唯一计算语义；matrix、operator、whitening、scaling、elimination 与 backend standard form 都是带 recovery map 的派生物。
- 后端能力族为 `SymmetricKktForm`、`ConvexQpForm` 与 `ConvexConicForm`；按数学能力路由，不按 Horizon、Normal 等领域名路由。
- Cubic 的 CPD reduction 必须证明 $\operatorname{range}(T)=\ker(P^\mathsf T)$ 且 $T^\mathsf TKT\succ0$；不能把原始 CPD Gram 当作 PSD QP Hessian。
- backend termination 只给候选结果；`RecoverAndVerify` 才能产生已求解模型。

参见 [canonical functional problem IR ADR](../adr/0001-canonical-functional-problem-ir.md)、[v1 kernel contracts ADR](../adr/0002-fix-v1-kernel-contracts.md) 和 [`CONTEXT.md`](../../CONTEXT.md)。

### 1.3 一手算法资料

- 对称不定分解的基线是 Bunch–Kaufman diagonal pivoting，而非无主元 LDL：[Bunch & Kaufman, 1977](https://doi.org/10.1090/S0025-5718-1977-0428694-0)；LAPACK 的 [`DSYTRF` source/spec](https://www.netlib.org/lapack/double/dsytrf.f)按该类方法生成 1×1/2×2 block diagonal factor。
- LAPACK [`DGEQP3`](https://www.netlib.org/lapack/double/dgeqp3.f)提供列主元 QR；[`DGESVD`](https://www.netlib.org/lapack/double/dgesvd.f)提供 singular values；[`DGELSD`](https://www.netlib.org/lapack/double/dgelsd.f)用 SVD 解最小二乘并按 `RCOND` 判 rank。
- Clarabel 的 homogeneous embedding 与原生 quadratic objective 见作者论文：[Goulart & Chen, 2024](https://arxiv.org/abs/2405.12762)。OSQP 的 ADMM、quasi-definite system 与 infeasibility detection 见作者论文：[Stellato et al.](https://arxiv.org/abs/1711.08013)。

算法资料定义能力，但不证明某个 Rust wrapper 已正确暴露全部量；wrapper 事实以下表及固定版本源码为准。

## 2. 稠密线性代数事实矩阵

| 候选（固定版本） | SPD / symmetric-indefinite | RRQR / SVD / condition | 依赖、平台、许可、MSRV | sparse / operator seam | 事实风险 |
|---|---|---|---|---|---|
| [`faer 0.24.4`](https://docs.rs/faer/0.24.4/faer/) | `Mat::llt`；`Mat::lblt` 为带 permutation 和 1×1/2×2 `B` 的 Bunch–Kaufman 型分解，factor 与 permutation 可读取 | `col_piv_qr` 暴露 `R`/`P`；full/thin SVD 与 ordered singular values；公开高层 API 未提供内建 numerical-rank 或 `rcond` 结论 | [README 与 manifest](https://docs.rs/crate/faer/0.24.4/source/)分别声明 pure Rust，以及 MIT/Rust 1.84；默认含 Rayon 与 sparse-linalg | 有 CSC sparse、sparse Cholesky/LU/QR，并有 [`matrix_free`](https://docs.rs/faer/0.24.4/faer/matrix_free/index.html) module | rank threshold、condition policy 与 inertia/near-zero-pivot acceptance 仍需项目实现；crate MSRV 较高 |
| [`nalgebra 0.35.0`](https://docs.rs/nalgebra/0.35.0/nalgebra/linalg/index.html) | `Cholesky`；[`LBLT`](https://docs.rs/nalgebra/0.35.0/nalgebra/linalg/struct.LBLT.html)明确实现 Bunch–Kaufman Algorithm A，`d()` 和 solve failure 可见 | `ColPivQR` 暴露 `R`/permutation，但 `is_invertible` 只查精确零；`SVD::rank(eps)` 由调用方提供阈值 | [manifest](https://docs.rs/crate/nalgebra/0.35.0/source/Cargo.toml)为 Apache-2.0、Rust 1.89；核心路径 pure Rust，支持 `no_std`/alloc 配置 | sparse feature 有 `CsMatrix`/`CsCholesky`，但 decomposition API 不是统一 operator 后端 | 当前 LBLT 是新近能力；需专门验证病态 KKT、inertia 和性能，不能从 API 丰富度推断中大型性能 |
| [`nalgebra-lapack 0.28.0`](https://docs.rs/nalgebra-lapack/0.28.0/nalgebra_lapack/) | LAPACK Cholesky；没有公开 symmetric-indefinite wrapper，因此需与 `nalgebra::LBLT` 或直接 LAPACK binding 组合 | [`ColPivQR`](https://docs.rs/nalgebra-lapack/0.28.0/nalgebra_lapack/colpiv_qr/struct.ColPivQR.html)调用 `?geqp3`，公开 rank 与 rank algorithm；另有 SVD | [manifest](https://docs.rs/crate/nalgebra-lapack/0.28.0/source/Cargo.toml)为 MIT，依赖 BLAS/LAPACK；未单列 rust-version，但依赖 Rust 1.89 的 nalgebra 0.35，故有效 floor 不低于 1.89 | dense LAPACK；无 sparse/operator abstraction | 默认 rank heuristic 是策略而非数学真理；native linking 与 binary distribution 成为应用责任 |
| [`ndarray-linalg 0.18.1`](https://docs.rs/ndarray-linalg/0.18.1/ndarray_linalg/) | `Cholesky`；[`SolveH`/`FactorizeH`](https://docs.rs/ndarray-linalg/0.18.1/ndarray_linalg/solveh/index.html)通过 `?sytrf`/`?hetrf` 做 Bunch–Kaufman | SVD/`?gesdd`；`LeastSquaresSvd` 通过 `?gelsd` 返回 singular values 与 rank；一般 LU 路径暴露 `?gecon` 的 `rcond` estimate；普通 `QR` wrapper 不做列主元 | [README](https://github.com/rust-ndarray/ndarray-linalg/blob/ndarray-linalg-v0.18.1/README.md)要求 OpenBLAS/Netlib/MKL 三选一并列出编译/系统依赖；MIT OR Apache-2.0；manifest 未声明 MSRV；README 的 tested matrix 仅承诺 x86_64 | dense ndarray + LAPACK；无 sparse/operator seam | Cargo feature additive，library 不应替最终应用锁 BLAS backend；`gelsd` 在 wrapper 中以 `rcond=-1` 使用 LAPACK 默认阈值，不能直接成为 GeoRBF rank policy |
| 直接 [`lapack 0.20`](https://docs.rs/lapack/0.20.0/lapack/) / `lapack-sys` | 可直接取 `?potrf/?sytrf/?sytrs/?sycon` 等缺失 routine | 可直接取 `?geqp3/?gesvd/?gesdd` 并自己保存 workspace/info | 最薄但 unsafe/低层；仍需选择 BLAS/LAPACK 实现并承担 ABI/整数宽度/线程配置 | dense Fortran ABI；operator seam 不由它提供 | 能力最全但 adapter、workspace、error mapping、layout、thread control 和安全测试工作最大 |

### 2.1 分解不等于诊断

**已验证事实**

- `faer` 和 `nalgebra` 都把 column-pivoted `R` 或 singular values 暴露给调用方；`nalgebra::SVD::rank(eps)` 明确把 `eps` 留给调用方。[faer decomposition docs](https://docs.rs/faer/0.24.4/faer/#matrix-decompositions)，[nalgebra SVD docs](https://docs.rs/nalgebra/0.35.0/nalgebra/linalg/struct.SVD.html)。
- `nalgebra-lapack` 的 RRQR rank 提供 fixed lower bound 与两个 scaled-epsilon heuristic；其源码本身警告 fixed bound 容易误判。[rank source](https://docs.rs/crate/nalgebra-lapack/0.28.0/source/src/colpiv_qr/rank.rs)。
- `ndarray-linalg` 的 `rcond` 是由 `?gecon` 估计的一范数 reciprocal condition number，而不是精确条件数；它适用于其一般 LU factorization path。[`ReciprocalConditionNum`](https://docs.rs/ndarray-linalg/0.18.1/ndarray_linalg/solve/trait.ReciprocalConditionNum.html)。

**推论**

- GeoRBF 应把 “decomposition computed” 与 “rank/condition accepted” 分开。rank tolerance 需要引用 scaled matrix、norm、dimension 与 expected mode；恢复诊断还要指回缺失的 physical/polynomial modes。
- Bunch–Kaufman factor 能支持 KKT solve 与 inertia/pivot inspection 的基础，但上述高层 crate 没有共同、版本稳定的 `rank_report()` 契约。若数值政策需要 inertia 或 reciprocal condition，可能要在 adapter 中分析 1×1/2×2 blocks、做 SVD oracle，或调用额外 LAPACK routine。

**未知/待测**

- 固定 GeoRBF matrix corpus 上的 backward error、pivot growth、false rank decisions 与线程确定性。
- `faer`/`nalgebra` pure-Rust kernels 与所选 OpenBLAS/MKL/Netlib 在目标 CPU/OS 的实际 crossover；上游 benchmark 不能替产品 workload benchmark。

## 3. Convex QP 与 SOCP 求解器事实矩阵

| 能力 | [`Clarabel.rs 0.11.1`](https://docs.rs/clarabel/0.11.1/clarabel/) | [`osqp.rs 1.0.1`](https://docs.rs/osqp/1.0.1/osqp/) |
|---|---|---|
| 数学范围 | 直接接受 PSD quadratic objective 与 zero/nonnegative/SOC 等 cone；同一路径覆盖 LP/QP/SOCP，[标准形](https://github.com/oxfordcontrol/Clarabel.rs/tree/v0.11.1#readme) | 只接受 $\tfrac12x^TPx+q^Tx,\ l\le Ax\le u,\ P\succeq0$；不支持 SOC，[Rust wrapper README](https://docs.rs/crate/osqp/1.0.1/source/README.md) |
| 输入/线性系统 | CSC；默认 Rust QDLDL，optional `faer-sparse`/PARDISO；当前 settings 要求 direct KKT solver | CSC；`osqp-sys` 静态编译 pinned OSQP 1.0.0 C source；safe settings 有 direct/indirect selector，但实际 compiled capability 仍须运行时/构建验证 |
| termination states | `Solved`, primal/dual infeasible, `Almost*`, max iterations/time, numerical error, insufficient progress, callback terminated，[enum](https://docs.rs/clarabel/0.11.1/clarabel/solver/enum.SolverStatus.html) | solved/inaccurate、primal/dual infeasible/inaccurate、max iterations/time、non-convex；Rust enum 没有 interrupt variant且 build script关闭 interrupt，[status](https://docs.rs/osqp/1.0.1/osqp/enum.Status.html)，[build script](https://docs.rs/crate/osqp-sys/1.0.1/source/build.rs) |
| residual / gap | public `solver.info` 有 primal/dual residual、primal/dual infeasibility residual、absolute/relative gap、objectives、κ/τ；solution 有 unscaled x/z/s，[DefaultInfo](https://docs.rs/clarabel/0.11.1/clarabel/solver/implementations/default/struct.DefaultInfo.html) | safe `Solution` 暴露 primal/dual residual、objective、x/y；[pinned OSQP C `OSQPInfo`](https://github.com/osqp/osqp/blob/236713ce9a56c182ac3230d52108f952afce1523/include/public/osqp_api_types.h#L128-L162)有 dual objective/gap，但 `osqp.rs 1.0.1` 没有对应 safe accessor，[Rust source](https://docs.rs/osqp/1.0.1/src/osqp/status.rs.html) |
| certificates | homogeneous embedding；infeasible statuses 明确说明返回 primal/dual certificate，`x/z/s` 可取 | `PrimalInfeasibilityCertificate::delta_y`、`DualInfeasibilityCertificate::delta_x` 是明确 safe types |
| scaling visibility | settings 可开关/限制 equilibration；public problem data 保存 `d/dinv/e/einv/c`，solution post-process unscale | settings 可控制 scaling iterations 与 scaled termination；safe wrapper 不暴露内部 scaling vectors；solution/certificates由 upstream unscale |
| gap/certificate acceptance | 可获得充分原始量，但 GeoRBF 仍须按 canonical units 复验；`Almost*` 不应静默并入成功/不可行 | certificate 可直接复验；缺 safe gap 是验收接口缺口，可由 GeoRBF 独立重算，或写受审 sys adapter |
| warm start / repeated solve | 支持固定 sparsity 的 P/q/A/b data update；但每次 `solve` 都调用 `default_start`，0.11.1 没有公开 initial iterate API，因此不可宣称 warm start，[solve source](https://docs.rs/crate/clarabel/0.11.1/source/src/solver/core/solver.rs) | 明确 `warm_start(x,y)`、x-only、y-only，并能在固定 sparsity 下更新 P/A/q/l/u；setting 可沿用上一解，[Problem API](https://docs.rs/osqp/1.0.1/osqp/struct.Problem.html) |
| condition / rank visibility | `LinearSolverInfo` 只给 solver name、threads、direct、KKT/factor nnz；没有 condition/rank/inertia | safe API 给 iteration/rho/timing/residual；没有 KKT condition/rank/inertia |
| regularization / refinement | settings 暴露 static/dynamic KKT regularization 与 iterative refinement；这是求解算法参数，但其影响仍须通过 canonical verification 隔离 | ADMM 的 `rho/sigma`、polishing regularization、iterative refinement 与 scaling 均可配置；不能解释为 GeoRBF FieldEnergy/ridge |
| 许可 / MSRV / native | [manifest](https://docs.rs/crate/clarabel/0.11.1/source/Cargo.toml.orig)：Apache-2.0、Rust 1.70；默认 QP/SOCP 路径无 BLAS/LAPACK C/Fortran，SDP/PARDISO features 会引入 native/licensing 变化 | [manifest](https://docs.rs/crate/osqp/1.0.1/source/Cargo.toml.orig)：Apache-2.0、Rust 1.63；`osqp-sys` 需要 CMake/C compiler 并静态构建 C library |
| 维护信号 | 0.11.1 发布线、2026 年仍有源码/CI 更新；MSRV workflow 验证 default features | 1.0.1 于 2025-04 发布并 pin OSQP C 1.0.0；Rust repo CI覆盖 Linux/macOS/Windows 与 Rust 1.63，但 wrapper 文档链接/changelog 仍残留 0.6.3 文本，需视为维护质量风险 |

### 3.1 certificate 与状态的解释

**已验证事实**

- Clarabel 的 `AlmostSolved` 与 `Almost*Infeasible` 使用一套 reduced-accuracy tolerance；`NumericalError` 与 `InsufficientProgress` 是独立 states，不是不可行证明。[settings](https://docs.rs/clarabel/0.11.1/clarabel/solver/implementations/default/struct.DefaultSettings.html)。
- OSQP 区分 exact/inaccurate infeasibility；upstream 文档给出 primal certificate $v$ 和 dual certificate $s$ 的代数验收条件。[OSQP infeasibility docs](https://osqp.org/docs/solver/index.html#infeasible-problems)，[status values](https://osqp.org/docs/interfaces/status_values.html)。

**推论**

- GeoRBF 不能仅把 solver enum 映射为领域 `Infeasible`。应先 unscale/recover certificate，再用派生标准形与 canonical provenance 独立检查有限性、归一化、Farkas/ray 条件及 tolerance；不能验证时只能保留 backend termination 或 numerical failure。
- dual infeasible 对 convex QP 通常对应 primal unbounded ray；GeoRBF 的产品诊断仍需区分 invalid convex reduction、unidentified field modes 与真正 unbounded standard form。
- Clarabel 能统一 QP/SOCP，不代表 equality-only KKT 应强制走 conic embedding；是否共享或分开 route 是 issue 8 的政策问题。

### 3.2 warm start 与 scaling

**已验证事实**

- Clarabel 可以复用 solver allocation并更新固定 sparsity data，但 0.11.1 的 `solve()` 在每次迭代前执行 `default_start()`；data update 不等于 iterate warm start。
- OSQP safe wrapper 提供完整/部分 warm start；其 scaling vectors 位于 C workspace 私有状态，safe wrapper 只暴露 controls 和 unscaled results。

**推论**

- 若 v1 problem snapshot 是一次性 batch solve，Clarabel 无 warm start 未必是 blocker；若未来做参数 sweep/交互更新，二者能力差异会成为实质 route 输入。
- 无论上游是否 unscale，GeoRBF 都应保存自己的 lowering/scaling/recovery map。solver 内部缩放只应当是额外的数值变换，不能成为物理语义或诊断来源。

## 4. native dependency、平台、许可、MSRV 与维护风险

| 维度 | 事实 | 决策输入 / 风险 |
|---|---|---|
| pure Rust floor | faer 1.84/MIT；nalgebra 1.89/Apache-2.0；Clarabel 1.70/Apache-2.0 | 选 pure Rust 仍会锁定 MSRV；nalgebra 0.35 的 Rust 1.89 对 library consumers 可能是明显成本 |
| LAPACK floor | nalgebra-lapack/ndarray-linalg 要求 BLAS/LAPACK provider；OpenBLAS/Netlib static 常需 C/Fortran toolchain，MKL 有额外许可条款 | 需决定 library 是否完全不选 provider、由最终 binary 选 feature；还需定义 oversubscription/thread determinism policy |
| OSQP native | `osqp-sys` 用 CMake + `cc` 静态构建 pinned C source，并按 target pointer width选择 index width | cross compile、WASM、musl、Android/iOS 等目标不能从桌面 CI 外推；必须在 GeoRBF target matrix 验证 |
| Clarabel optional native | default QP/SOCP 不启用 SDP/BLAS/PARDISO；相关 features 会改变依赖和许可 | v1 只需 SOC 时应避免无意启用 SDP/PARDISO；Cargo feature unification 要纳入 audit |
| maintenance | 四个主候选均有版本化源码；faer/nalgebra/ndarray-linalg/Clarabel 在 2026 有更新，osqp.rs 最近 release 为 2025-04 | 活跃并不等于 API 稳定；应 pin minor/patch、保存 backend fingerprint，并用 acceptance corpus 管 upgrade |
| MSRV unknown | ndarray-linalg 与 nalgebra-lapack package manifest 未显式声明 rust-version；后者的 nalgebra dependency给出有效下界，前者没有可审计声明 | “能在当前 CI 编译”不是 MSRV contract；若准入必须由 GeoRBF CI 自己固定最低 toolchain |

### 4.1 简短边界候选

- **POUNCE 0.6（watchlist）**：**已验证事实**：[v0.6.0 于 2026-06-20 发布](https://github.com/jkitchin/pounce/releases/tag/v0.6.0)，EPL-2.0、默认 pure Rust；[project README](https://github.com/jkitchin/pounce/blob/v0.6.0/README.md) 声明 SOCP/certificate/warm start，而 tagged `pounce-convex` [source](https://github.com/jkitchin/pounce/blob/v0.6.0/crates/pounce-convex/src/ipm.rs) 已公开相关 entries，但 [crate README](https://github.com/jkitchin/pounce/blob/v0.6.0/crates/pounce-convex/README.md) 仍称 LP/QP 之外的 conic family 在建设。**推论**：这是很新的 LP/QP/SOCP challenger，而非成熟度证据。**未知**：稳定 diagnostics/API、MSRV 与平台 contract。
- **Totsu**：**已验证事实**：[`totsu_core`](https://docs.rs/totsu_core/0.1.1/totsu_core/) 是 operator-oriented first-order conic solver，含 pure-Rust 路径与 SOC；[upstream](https://github.com/convexbrain/Totsu) 于 2026-03-07 archive。**推论**：archive 是实质 maintenance 边界。**未知**：后续维护与兼容性 contract。
- **MOSEK Rust API**：**已验证事实**：官方 API 覆盖 [QP](https://docs.mosek.com/latest/rustapi/tutorial-qo-shared.html) 与 [SOCP](https://docs.mosek.com/latest/rustapi/tutorial-cqo-shared.html)，并要求 [native installation](https://docs.mosek.com/latest/rustapi/install-interface.html)；Optimizer 是需有效 license 的 [commercial product](https://docs.mosek.com/latest/licensing/intro.html)。**推论**：它是 proprietary/native/commercial comparator，不能假设为可自由再分发的默认 backend。**未知**：GeoRBF 分发与部署所需授权。
- **HiGHS**：官方能力是 LP/MIP/convex QP，不是 SOCP；Rust binding 还会引入 C++ native library。因此它最多是另一条 QP route，不填补 DirectedNormal SOCP，[HiGHS project](https://github.com/ERGO-Code/HiGHS)，[`highs` crate](https://docs.rs/highs/latest/highs/)。本次未深挖其 residual/certificate safe API。
- **ECOS**：上游是 C SOCP solver，但 crates.io 的 [`ecos 0.1.0`](https://docs.rs/ecos-rs/latest/ecos/) 基本是裸 FFI surface，缺本票所需的安全 status/residual/certificate adapter 证据；不作为主候选。其后续维护与平台状态本次未知。
- **SCS**：上游支持 cone problems与 certificates，但官方接口列表没有 Rust；未找到与 Clarabel/OSQP 同级、可审计的官方 Rust package，[SCS API](https://www.cvxgrp.org/scs/api/index.html)。若未来引入第三方 wrapper，应重新做版本、FFI ownership、status mapping 与 certificate audit。

这些行是调查边界，不是永久拒绝决定。

## 5. 规模、内存与复杂度

令：

- $n$：backend primal variables，包括 representer coordinates、polynomial/semantic latents 与 slacks；
- $m$：affine/cone rows；
- $k$：某次直接法实际 factor 的 KKT order；
- `nnz(data)` 与 `nnz(L)`：输入 CSC 和 factor after-fill nonzeros。

### 5.1 稠密线性场景

**算术事实**

- 一个 full dense $N\times N$ f64 buffer 需 $8N^2$ bytes；$N=2{,}000$ 时约 30.5 MiB。
- 只存一侧三角的 values 是 $8N(N+1)/2$ bytes，$N=2{,}000$ 时约 15.3 MiB；但本表高层 dense matrix types通常持有 full rectangular buffer，factor、RHS、copy 和 workspace 会形成多个 buffer。
- dense Cholesky/LDL/QR/SVD 都有 $O(N^3)$ factorization time；Cholesky 主项约 $N^3/3$，$N=2{,}000$ 约 2.67×10⁹ floating operations。列主元 QR/SVD 常数更大，适合作 rank oracle 或 preflight，不应假设与一次 Cholesky 同价。

**推论**

- $N=2{,}000$ 的一个 buffer 并不危险，但装配矩阵 + factor copy + RRQR/SVD workspace + recovery artifacts 很容易进入数百 MiB。后续性能门槛必须测 peak RSS，不只测输入矩阵大小。
- 默认 Cubic 的 augmented KKT 可能大于 constraint count；CPD null-space reduction 也可能生成 dense $T$ 或 dense reduced Hessian。显式 $Z$ 的内存/乘法成本是 issue 8 的算法决策，不能由 crate 表替代。

### 5.2 QP/SOCP 场景

**算术事实**

- 若 dense PSD $P$ 以 CSC upper triangle存储，value count 为 $n(n+1)/2$。仅按 64-bit value + 64-bit row index 粗算，$n=2{,}000$ 约 30.5 MiB，尚未计 A、column pointers、KKT assembly 与 fill。
- sparse direct solver 的实际内存更接近 `nnz(data) + nnz(L)`；Clarabel 的 `LinearSolverInfo` 正好暴露 KKT 与 factor nnz，但不暴露 condition/rank。
- interior-point 每次迭代需要 KKT update/factor/solve；ADMM/OSQP setup factorization 后进行多次较便宜 iteration，并在 rho/data matrix update 时可能 refactor。二者不能只用 iteration count横比。

**推论**

- GeoRBF 的 kernel Gram/Hessian 很可能 dense，即使 A 与 SOC blocks sparse，也会削弱 generic sparse-solver 的内存优势；Wendland C2 才可能形成真正 sparse pairing，但具体 sparsity 依 support radius 与 functional coupling。
- “最多 2,000 lowering constraints” 只能是 workload selector。每次 solve 应记录
  $(n,m,k,\operatorname{nnz}P,\operatorname{nnz}A,\operatorname{nnz}L,\text{cone blocks})$，再决定是否允许 dense、sparse direct 或未来 operator route。

**未知/待测**

- GeoRBF 代表性 equality/QP/SOCP corpus 上的 KKT order、fill ratio、peak memory、setup/solve time 与 accuracy distribution。
- 多线程 BLAS/faer 对 2,000 级问题是否抵消调度与 oversubscription；需要统一线程数再 benchmark。

## 6. 给后续数值政策票的决策输入

本节只列需要决定的轴，不给最终选择：

1. **Dense linear baseline**：pure Rust Bunch–Kaufman/Cholesky 是否足够，还是要用 LAPACK 作为 production route 或 oracle；若组合多个 crate，矩阵转换和重复分配预算是多少。
2. **Rank/condition policy**：RRQR 还是 SVD 为权威 oracle；tolerance 如何随 scaling、matrix norm、dimension、expected CPD modes变化；是否需要 inertia/`rcond` 及其 backend adapter。
3. **QP/SOCP route**：QP 是否与 SOCP 都走 Clarabel，还是 QP 另走 OSQP；route 的理由应是 certificate/warm-start/accuracy/workload，而不是领域类型。
4. **Status acceptance**：`Almost*` / inaccurate 是否只保留为 backend diagnostics；certificate 独立验收条件；max iteration/time/insufficient progress/numerical error 的精确映射。
5. **Scaling 与 regularization**：GeoRBF 自有 scaling、solver equilibration、KKT static/dynamic regularization、OSQP rho/sigma/polish 各自 provenance；禁止把算法 regularization解释成 FieldEnergy 或 hidden ridge。
6. **Warm start**：一次性 v1 是否明确不承诺 warm start；若承诺，snapshot compatibility、fixed sparsity、candidate verification 与 deterministic fallback 怎么定义。
7. **Dependency budget**：MSRV、native toolchain、BLAS provider、threading、cross-compile target、许可与 binary size 的可接受边界。
8. **Future seam**：Canonical IR/AlgebraicPlan 保持 storage independent；v1 backend adapter 可以 dense/CSC specific，但不能把 faer/nalgebra/ndarray/Clarabel/OSQP types 写入 canonical semantics。

## 7. 研究限制与复现实验清单

本研究是固定版本的静态 capability audit，没有编译或 benchmark 各候选。以下事项必须在准入前用 GeoRBF 自有 fixture 验证：

- Linux/macOS/Windows 目标与最低 Rust toolchain 的 clean builds；若用 LAPACK，再覆盖每个承诺 provider；
- SPD、near-SPD、singular、symmetric-indefinite、scaled KKT、known-rank CPD polynomial blocks；
- feasible、primal infeasible、dual infeasible/unbounded、nearly infeasible、badly scaled、max-iteration 与 numerical-failure QP/SOCP；
- solution、dual、gap 和 certificate 经 GeoRBF recovery 后的 independent recomputation；
- $n/m/k$ 分层到 2,000-constraint baseline 的 time、peak RSS、factor fill 与 repeatability；
- dependency license、Cargo feature unification、unsafe/FFI surface 与 supply-chain audit。

未验证事实不得升级为产品承诺：尤其是平台覆盖、rank threshold 可靠性、solver status 到领域诊断的直接映射，以及任一候选在 GeoRBF dense kernel workload 上的相对性能。
