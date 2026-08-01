# Cubic CPD KKT/QP/SOCP recovery spike

This unpublished, throwaway crate is the T06/T07 experiment for [GeoRBF issue #14](https://github.com/qingsongyukuai/GeoRBF/issues/14). It is outside the future product crate tree and does not admit any production dependency.

The manufactured canonical case uses ten single-support generalized functionals

```text
F(f) = a f(x) + b^T grad f(x)
```

including ordinary values, three coordinate derivatives, a value/derivative contraction, and a directional derivative. It assembles the physical Cubic pairing `K`, constructs the complete `Pi_1 = span{1,x,y,z}` pairing in the prescribed normalized coordinates, and retains an explicit recovery map to physical polynomial coefficients.

The experiment runs one manufactured physical truth through three solver forms:

- Equality: a symmetric augmented KKT containing the full `K`, full `P`, CPD side conditions, shared-level latent/gauge relations, and all hard equalities;
- QP: the same canonical equalities plus an inactive affine bound, reduced through an implicit Householder QR null-space operator;
- SOCP: the same canonical equalities plus an interior second-order cone, using the same reduction and recovery.

All standard forms receive eight rounds of max-norm Ruiz scaling. Factors are nearest powers of two, each round is clipped to `[2^-8,2^8]`, cumulative factors to `[2^-32,2^32]`, KKT scaling is diagonal congruence, and the SOC block receives one common row factor. Clarabel's own equilibration remains enabled and is an additional recorded backend implementation detail, not physical semantics.

Replay with the exact Rust 1.85 toolchain and locked temporary backends:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo run --locked
```

The integration tests exercise only the public evidence seam. They check complete polynomial rank, the CPD side condition, reduced strict positivity and symmetry, affine reproduction, KKT inertia/backward error, independent QP/SOCP residuals, physical recovery, cross-route canonical agreement, and the three required counterexamples. No route deletes a polynomial mode, adds ridge/jitter, anchors a field automatically, or changes kernel.
