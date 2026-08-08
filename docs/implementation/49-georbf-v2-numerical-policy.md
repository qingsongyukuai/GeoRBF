# Issue 49: unique georbf-v2 numerical policy

Issue: [#49](https://github.com/qingsongyukuai/GeoRBF/issues/49)

Primary evidence seam: the public `FitConfiguration` -> `ProblemSnapshot` ->
`fit` -> `FitReport` and `RepresentationEvidence` workflow.

## One executable policy

`FitConfiguration::default` resolves directly to `georbf-v2`, and every
snapshot retains that complete configuration. Callers can inspect the policy
identity and change only the thread budget: there is no runtime policy
selector, fallback, or compatibility switch. `NumericalPolicyId::georbf_v1`
remains readable as the immutable historical v0.2.0 identity, but cannot be
installed into a fit configuration.

The v2 constant is the single internal owner of scaling, decomposition and
certificate limits, precision rescue, canonical recovery, and reliable query
acceptance. Equality KKT and Convex QP reports both carry that same identity;
backend fingerprints and factorization settings remain separate evidence.

## Compatibility and audit evidence

The v0.2.0 construction, relation, fit, and query signatures remain unchanged.
Hardness and soft-loss semantics still enter the canonical solver form from the
original domain inputs, complete Pi1 remains the Cubic polynomial contract, and
physical recovery continues to use
`1e-10 * characteristic_scale + 1e-8 * relation_reference_scale` per relation.

`RepresentationEvidence` now reports the hard-to-soft conversion count
explicitly. The only v2 route reports zero, alongside the existing zero
truncated-mode and no-problem-regularization evidence. Backend factorization
regularization remains recorded per attempt and never replaces verification of
the unregularized canonical problem.

The public regression suite covers the default and snapshot identities,
historical v1 readability, the rescued small-positive-mode case, successful
canonical recovery, invalid input, interpretable rank deficiency,
infeasibility, backend-contract failure, and recovery-verification failure.
These checks compare domain diagnoses and canonical observables, not bitwise
coefficients, FieldEnergy, or query samples.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --test public_problem_contract
cargo test --locked --test public_fit_diagnostics
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
```
