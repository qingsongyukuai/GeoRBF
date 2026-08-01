# Cubic CPD KKT/QP/SOCP recovery spike

This unpublished, throwaway crate is the T06/T07 experiment for [GeoRBF issue #14](https://github.com/qingsongyukuai/GeoRBF/issues/14). It is outside the future product crate tree and does not admit any production dependency.

The manufactured canonical case uses ten single-support generalized functionals

```text
F(f) = a f(x) + b^T grad f(x)
```

including ordinary values, three coordinate derivatives, a value/derivative contraction, and a directional derivative. It assembles the physical Cubic pairing `K`, constructs the complete `Pi_1 = span{1,x,y,z}` pairing in the prescribed normalized coordinates, and retains an explicit recovery map to physical polynomial coefficients.

The experiment defines one physical canonical relation set: eleven equalities, one affine upper bound, and one SOC. Because the equalities fix every functional value and the semantic latent, the bound and cone are provably redundant. Three algebraic plans preserve the same feasible set while retaining the constraints needed to trigger each solver family:

- Equality: elides both proven-redundant convex relations and uses a symmetric augmented KKT containing the full `K`, full `P`, CPD side conditions, shared-level latent/gauge relations, and all hard equalities;
- QP: retains the affine bound and elides only the proven-redundant cone, reduced through an implicit Householder QR null-space operator;
- SOCP: retains the cone and elides only the proven-redundant affine bound, using the same reduction and recovery.

The owned `risk-spike-14-v1` numerical policy drives rank, convexity, inertia, recovery, scaling, and acceptance thresholds. All standard forms receive eight rounds of max-norm Ruiz scaling. Factors are nearest powers of two, each round is clipped to `[2^-8,2^8]`, cumulative factors to `[2^-32,2^32]`, KKT scaling is diagonal congruence, and the SOC block receives one validated common row factor. Clarabel's own equilibration remains enabled and is an additional recorded backend implementation detail, not physical semantics.

Replay with the exact Rust 1.85 toolchain and locked temporary backends:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo run --locked
```

The integration tests exercise only the public evidence seam. They check complete polynomial rank with acceptance/rejection/gray-zone bands, the CPD side condition, reduced strict positivity and symmetry, affine reproduction, KKT inertia/backward error, independent QP/SOCP residuals (including zero-cone equality slack), physical slack equations, physical recovery, cross-route canonical agreement, and the three required counterexamples. No route deletes a polynomial mode, adds ridge/jitter, anchors a field automatically, or changes kernel.
