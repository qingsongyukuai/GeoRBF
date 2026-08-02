# Issue 18: Cubic Equality NumericalPolicy and Recover and Verify

Issue: [#18](https://github.com/qingsongyukuai/GeoRBF/issues/18)

Evidence seams: T06, T07, T08, T12

Requirements: IR-001, IR-007, IR-009, IR-014; NUM-003, NUM-004,
NUM-006–NUM-012, NUM-014, NUM-015; DIA-003, DIA-005; VAL-005

## Physical canonical form and reversible solve coordinates

The fitting uses and provenance remain in the physical Canonical Problem IR.
The derived Cubic form uses the deterministic similarity transform
`x_hat = (x - c) / L`, where `c` is the axis-aligned bounding-box center and
`L` is the maximum Euclidean support radius. A zero extent records explicit
degenerate-extent evidence and resolves `L` to one. A non-finite length or a
Cubic `L^3` recovery scale that cannot be represented in `f64` is rejected
before assembly or backend invocation.

The derived functional keeps value coefficients and divides physical-gradient
coefficients by `L`. Cubic pairing and the complete `Pi1` polynomial pairing
are assembled only in these solve coordinates. Recovery applies the exact
inverse maps

```text
c_physical = c_hat / L^3
b_physical = [b0 - (b_linear / L) dot c, b_linear / L]
```

and independently maps both field and polynomial coefficients back to solve
coordinates with a `1e-11` round-trip limit. Backend residual/tolerance vectors
have an explicit forward/inverse map through the algebraic scaling factors.

## Rank, scaling, inertia, and attempts

Exact zero structure is rejected before numerical scaling. General full-rank
decisions retain a column-pivoted RRQR screening ratio and confirming SVD
spectrum. For `d = max(m, n)` and `rho = sigma_min / sigma_max`, the owned
`georbf-v1` policy rejects at `rho <= 64 u d`, accepts at
`rho >= 4096 u d`, and returns typed Numerical Decision Gray Zone evidence
between the bands. Cubic polynomial rank applies the same structure-first,
RRQR-then-SVD policy before any field solve.

Every valid symmetric KKT receives exactly eight Ruiz max-norm rounds. Each
diagonal-congruence factor is quantized to the nearest `2^k`, clipped to
`[-8, 8]` per round and `[-32, 32]` cumulatively. The report retains all round
exponents, cumulative inverse factors, saturation evidence, and exact
matrix/RHS/residual/tolerance recovery maps. Symmetry is preserved by applying
the same factor to each row and corresponding column.

The scaled KKT must be confirmed full rank and have expected convex Equality
inertia `(primal, equality, zero) = (n, m, 0)` before a backend candidate can
exist. The deterministic Attempt Plan is:

1. faer LBLT Bunch-Kaufman with at most two refinement corrections;
2. one faer SVD rescue, only after the SVD rank evidence confirms full rank.

Each executed attempt records its sequence, algorithm and resolved faer
fingerprint, requested/actual threads, scaling summary, refinement count,
termination, normalized backward error, certificate absence, and failure
reason. Candidate finiteness and scaled normalized backward error `<= 1e-11`
use the same contract on both attempts. NumericalPolicy identity remains
separate from BackendFingerprint.

## Recover and Verify

An accepted backend-standard-form candidate is not a model. The independent
recovery boundary verifies the inverse coordinate map and every usage-edge
provenance entry before reconstructing the physical Cubic field. It then
recomputes, without reusing the backend equality rows:

- physical field and full gradient observations;
- every dimensioned hard-equality residual;
- the complete Cubic side condition;
- polynomial and field coefficient round trips;
- native Cubic FieldEnergy and its solve-coordinate round trip;
- the all-hard total objective `0.5 * FieldEnergy`;
- finiteness of all recovered quantities.

This tracer has no semantic latent input, so the verified recovered latent
count is zero. Backend-standard-form rejection and Recovery Verification
Failure are different typed boundaries, and neither failure contains a
candidate, partial field, or best-effort model. Corrupted coordinate and
provenance recovery maps have dedicated fail-closed corpus cases.

## Checked peak plan

The issue #16 capacity plan is extended before any allocation or backend call
to cover RRQR auxiliaries, singular-value and inertia analysis, full SVD rescue
`U/V` storage and exact faer workspaces, recovery buffers, and the expanded
attempt/scaling report. The deterministic all-live conservative boundary now
accepts KKT dimension 9,841 and rejects adjacent dimension 9,842 at
8,590,810,592 planned bytes, 876,000 bytes over the 8 GiB limit. This replaces
the smaller issue #16 evidence envelope because the bounded rescue path is now
part of the resolved policy.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```

The accumulated suite includes manufactured physical recovery, rank-band gray
zones, bad polynomial rank, wrong inertia, non-finite candidates, bounded SVD
rescue exhaustion, Ruiz reversibility, damaged recovery maps, provenance
round-trip failure, and adjacent 8 GiB capacity evidence.
