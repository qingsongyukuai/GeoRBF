# Issue 23: Equality structured diagnostics and Fit Report

Issue: [#23](https://github.com/qingsongyukuai/GeoRBF/issues/23)

Evidence seams: T08, T12

Requirements: PAPI-013, PAPI-017; DOM-003, DOM-019, DOM-020; IR-009,
IR-014; DIA-001–DIA-003, DIA-005–DIA-009; NUM-007–NUM-012, NUM-014

## Strict fit outcome and complete report

`ProblemSnapshot::fit` remains a synchronous
`Result<FitSuccess, FitFailure>`. A success owns exactly one immutable accepted
model and its `FitReport`. A failure owns a `ProblemDiagnosis` and the same
report type; neither its type nor its accessors contain a field, candidate, or
best-effort model. Build errors remain `BuildFailure` values from the mutable
builder phase and never masquerade as fit diagnoses.

The report retains the evidence available at every reached boundary: resolved
problem contract and dimensions, canonical hard-relation assessments, shared
level values, Cubic polynomial/reduced-pairing analysis, rank and inertia,
bounded attempts and backend fingerprints, physical canonical acceptance,
capacity, preflight proof, and recovery rejection evidence. Success continues
to report FieldEnergy, total objective, shared values, and every original-unit
residual. Exact hard duplicates share a canonical equality while preserving a
separate assessment for each caller `SourceId`; no public per-source dual is
invented.

## Aggregate preflight evidence and primary diagnosis

Preflight accumulates all safely obtainable proof before choosing the primary
diagnosis. It no longer returns after the first disconnected singleton or
first exact conflict. Every `UninformativeSharedLevelSetEvidence`, direct exact
conflict, and relation-graph conflict is retained. Conflict lists sort by
semantic role and then caller identity; singleton evidence sorts by `GroupId`
and member `SourceId`. Singular accessors remain shorthand for the first item,
while plural accessors expose the complete ordered corpus.

The v0.1.0 Equality priority is explicit rather than implied by control flow:

1. unresolved semantics (`UninformativeSharedLevelSet`);
2. `DirectInputConflict`;
3. unidentified additive or field mode;
4. `CapacityExceeded`;
5. `BackendContractViolation`;
6. `RecoveryVerificationFailure`;
7. unclassified numerical or limit failure.

Thus a shift-invariant contradictory gradient problem reports the direct
conflict while retaining gauge evidence. A shift-invariant over-capacity
problem reports the gauge problem while retaining the checked 8 GiB capacity
plan. No large allocation or backend call occurs in either case.

Capacity planning has two checked stages. A source-lifecycle guard is computed
before evidence assembly and charges 2 KiB of fixed lowering/report storage per
scalar relation plus eight copies of caller identifier bytes. Lowering and the
eventual success report are separate 1 KiB/four-copy phases in that sum. When
the combined lifecycle is over budget but the lowering-only subplan is within
budget, diagnosis is deferred until that linear canonical/preflight pass has
recovered any higher-priority conflict or field-mode proof. If even the
lowering-only subplan exceeds 8 GiB, the fit returns without attempting it. No
dense KKT allocation, factorization workspace, or backend call is attempted in
either case. After deduplication, the all-live plan uses distinct solver-row,
canonical-relation, and source-report sizes, so exact duplicates consume linear
audit storage but do not inflate dense KKT storage. Once lowering has
completed, every preflight failure reports the exact canonical and KKT
dimensions. If an exact capacity failure coexists with a recoverable `Pi1` null
mode, the bounded four-column polynomial analysis retains both proofs and
selects the higher-priority `UnidentifiedFieldMode` diagnosis. If that bounded
analysis reaches a gray zone or another typed analysis failure, the report
retains it as secondary evidence while capacity remains primary.

## Interpretable rank deficiency

Low rank has semantic meaning only after its null mode is recovered to a
canonical field concept. Cubic polynomial failure now handles rectangular
pairings explicitly, obtains a deterministic right singular vector when an
exact missing column is not enough, maps the mode through the solve-coordinate
recovery, and independently checks its normalized canonical residual and
round trip. A reduced-pairing rank loss is evidence about a derived numerical
matrix only; until a physical field mode is independently constructed and
evaluated through Canonical IR, it remains `NumericalFailure` rather than being
given field semantics.

Only a mode passing the unchanged `georbf-v1` recovery limit produces
`InterpretableRankDeficiencyEvidence` and an `UnidentifiedFieldMode` diagnosis.
The public proof names the canonical Cubic concept, algebraic domain, stable
source/role provenance, verified residual, backend-invocation state, and the
fact that no hidden regularization was applied. A generic KKT rank loss without
such a canonical recovery remains `NumericalFailure`; it is not guessed to be
an unidentified field.

## Termination, backend contract, and recovery

`SolveAttemptTermination` now records only how an attempt stopped. The faer
Equality adapter reports `CandidateProduced` whether the candidate later passes
or fails its backend-standard-form checks; structured attempt failure evidence
records non-finiteness or excess backward error separately. Future reduced
accuracy, infeasibility-candidate, limit, insufficient-progress, callback, and
numerical terminations have distinct public variants and do not imply a
`ProblemDiagnosis`.

Backend-standard-form rejection remains `BackendContractViolation`. A candidate
that passes that adapter boundary but fails canonical recovery, round trip,
provenance, side condition, or physical relation acceptance remains
`RecoveryVerificationFailure`. Both paths retain their attempts and evidence,
and neither returns a model. No ridge, jitter, polynomial deletion, automatic
anchor, tolerance relaxation, or kernel fallback is introduced.

Reliable general Farkas and recession certificates remain outside this
Equality milestone, so no unsupported infeasible or unbounded conclusion is
fabricated.

## Validation evidence

`tests/public_fit_diagnostics.rs` exercises T08 and T12 through the public fit
boundary. It verifies permutation-stable aggregate preflight evidence,
diagnosis priority with retained secondary evidence, checked capacity before
allocation, a nontrivial rectangular `Pi1` unidentified mode, and the separation
between candidate termination and canonical acceptance. Test-only backend
boundary injection runs through `ProblemSnapshot::fit` and verifies that a
rejected candidate maps to backend-contract diagnosis while retaining
`CandidateProduced`; an uninterpreted KKT rank loss remains numerical.

The accumulated public/property corpus continues to cover exact duplicate
multi-source retention, graph conflicts, additive gauge invariance, recovery
failure, success report observables, and input permutations. Existing internal
failure injection continues to cover rank gray zones, negative curvature,
capacity, backward-error and non-finite candidate rejection, scaling damage,
provenance damage, and physical recovery damage without freezing basis, row,
pivot, coefficient, or iteration layouts.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --test public_fit_diagnostics
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```
