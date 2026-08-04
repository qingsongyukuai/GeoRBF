# Issue 35: Convex Relations verification and diagnosis closure

Issue: [#35](https://github.com/qingsongyukuai/GeoRBF/issues/35)

Primary evidence seams: T01, T04–T08, T10–T12, T14–T16

Requirements: PAPI-013, PAPI-017; IR-008–IR-009, IR-014;
NUM-007–NUM-015; DIA-001–DIA-009; VAL-005, VAL-007–VAL-008,
VAL-015

## One candidate and evidence policy

Every public v0.2 Convex Relation lowers to the single physical
`CubicCanonicalProblem` and uses the capability-selected QP route. Clarabel
`Solved` and `AlmostSolved` terminations are candidate-producing events only.
Both pass the same independently recomputed scaled primal, dual-cone,
stationarity, complementarity, and relative-gap limit of `1e-8`. Recovery then
maps primal, slack, and dual values through the exact GeoRBF scaling and
recomputes the same five-part envelope in the physical solver-independent QP
form before the unchanged canonical Recover-and-Verify boundary. A successful
`AlmostSolved` attempt remains visible as
`ReducedAccuracyCandidateProduced`; no tolerance is widened.

Limit, insufficient-progress, callback, numerical, and unverified certificate
terminations never contain a public model or canonical acceptance evidence.
The fixed attempt plan executes Standard and at most one Robust profile. Both
profiles reuse the same solver-independent form, objective, relation hardness,
kernel, loss semantics, numerical policy, and canonical tolerance. Only the
fully fingerprinted backend numerical settings change.

If multiple executed attempts ever carry independently validated but different
candidate, infeasible, or unbounded conclusions, the executor rejects the
sequence as `NumericalConsistencyFailure`. It never selects a conclusion by
attempt order. This defensive rule is separate from an unverified status retry.

## Farkas and recession verification

A primal-infeasibility termination becomes `InfeasibleProblem` only after the
Clarabel dual ray is recovered through the exact GeoRBF row scaling. The
solver-independent verifier normalizes it to unit infinity norm and checks all
quantities are finite, `A^T z`, the dual cone, and strict normalized
`-b^T z`. Residual and cone violations must not exceed `1e-8`; separation must
be at least `1e-7`; the scaled-ray round trip must not exceed `1e-11`; and the
complete canonical-to-form provenance map must have independently passed.

A dual-infeasibility termination follows the symmetric policy. The Clarabel
primal ray is recovered through variable scaling and normalized. Independent
verification requires zero quadratic curvature `P d`, equality recession,
inequality-cone feasibility, and strict objective descent `-q^T d`, with the
same `1e-8`, `1e-7`, and `1e-11` limits. Only then is
`ProblemDiagnosis::UnboundedProblem` formed. Invalid or absent Farkas and
recession rays remain numerical failures after the bounded retry plan.

The public `InfeasibilityCertificateEvidence` and `RecessionRayEvidence`
expose normalized checks, fixed limits, recovery error, provenance status, and
stable source/group/semantic-role associations. They never expose raw rays,
Clarabel enums, rows, columns, cones, or non-unique dual values.

## Complete canonical acceptance and recovery evidence

Successful QP reports expose both the independently recomputed scaled backend
envelope and the five-part residual envelope recomputed after physical
canonical recovery. `CanonicalAcceptanceEvidence` also contains:

- backend-standard-form verification and residual;
- maximum physical hard affine-inequality violation;
- complete physical Cubic side condition and hard equality maxima;
- GeoRBF scaling, Householder reduction, polynomial, field coefficient,
  FieldEnergy, whitening, objective, tolerance, and Clarabel internal-scaling
  round trips;
- provenance and finiteness decisions.

All recovery/scaling round trips retain the fixed `1e-11` acceptance boundary.
`CubicAnalysisEvidence` now reports both ends of the accepted reduced spectrum
and a condition estimate, complementing the existing complete polynomial rank
evidence without exposing a basis or pivot layout.

QP recovery failures publish the reached backend-standard-form residual,
reduction and scaling errors, deterministic rejection reasons, and stable
canonical source associations. They explicitly report that no model was
produced. FieldEnergy, total objective, shared values, relation residuals,
slacks, violations, active state, and group-only covariance contribution remain
available only after complete canonical acceptance. Exact hard duplicates keep
all caller assessments but do not invent per-source duals; covariance members
do not receive invented per-member objective contributions.

## Diagnosis precedence

The completed Convex priority is:

1. unresolved semantics, including unresolved Axial polarity;
2. direct local, geometric, or relation-graph conflict;
3. unidentified additive or observable field mode and interpretable rank loss;
4. checked capacity rejection;
5. independently validated infeasible or unbounded proof;
6. backend-standard-form contract violation;
7. canonical recovery verification failure;
8. numerical consistency, unverified termination, limit, or other numerical
   failure.

Preflight returns before backend execution when it owns a higher-priority
proof. Later execution diagnoses retain every attempt and all safely available
secondary evidence. No failure type has a field, candidate, coefficient,
partial model, or best-effort query surface.

## Cumulative corpus

`tests/public_convex_relations_closure.rs` constructs Field Value Bounds,
Directional Derivative Intervals, Younger Than, Field Level Order, Field
Separation Interval, Point-to-Level-Set, Directed Normal, explicitly resolved
Axial Normal, and a covariance group in one manufactured affine problem. It
checks their deterministically ordered physical reports, shared values,
group-only covariance objective, FieldEnergy, total objective, complete QP
residuals, condition evidence, canonical acceptance, and query-capable model.

The accumulated public suites supply relation-specific success, soft loss,
duplicate, permutation, frame/unit covariance, direct conflict, and validated
general infeasibility cases. Contract tests add accepted Almost Solved,
unverified limit exhaustion, backend-contract corruption, scaling/reduction/
objective/provenance recovery corruption, invalid Farkas and recession rays,
validated unboundedness on a manufactured solver-independent QP, and
contradictory-conclusion rejection. Regressions assert only canonical
observables and stable provenance, never coefficients, basis, row layout,
pivots, iterations, raw duals, or backend trajectory.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --test public_convex_relations_closure
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```
