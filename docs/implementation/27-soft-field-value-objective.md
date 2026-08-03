# Issue 27: Soft Field Value and physical objective

Issue: [#27](https://github.com/qingsongyukuai/GeoRBF/issues/27)

Primary evidence seams: T01, T03, T04, T07–T11

Requirements: PAPI-010–PAPI-013; DOM-017–DOM-018; IR-003, IR-005–IR-008,
IR-014; KER-008; NUM-005, NUM-009–NUM-010; DIA-001, DIA-006, DIA-008–DIA-009;
VAL-004–VAL-005, VAL-008, VAL-015

## Public boundary

`FieldValueObservation::try_new` remains the hard constructor. The two named
soft constructors accept distinct crate-owned checked types:

- `try_with_quadratic_penalty(..., QuadraticPenalty)` records a finite positive
  optimization weight with no statistical interpretation.
- `try_with_standard_deviation(..., StandardDeviation)` records finite positive
  statistical uncertainty and derives its canonical precision as `1 / sigma^2`.

There is no public universal enforcement or loss enum. The successful report
retains the distinction through the typed `quadratic_penalty` and
`standard_deviation` accessors on each `SoftFieldValueAssessment`.

`FieldEnergyNormalization::try_new` accepts only finite positive factors.
`ProblemBuilder::set_field_energy_normalization` stores the single explicit
normalization. An all-hard build still resolves an omitted value to exactly
one; any soft Field Value without it contributes
`BuildError::MissingFieldEnergyNormalization` to the deterministic aggregate
and `BuildFailure::into_builder` retains the repairable builder.

## Canonical and numerical path

Lowering keeps hard equalities, soft equality residual channels, and soft
objective terms as separate canonical collections. Hard duplicates may still
share one solver equality while retaining every source assessment. Soft
duplicates are never merged: every SourceId owns a residual, typed loss, and
objective contribution even when structural functionals are interned in the
representer span.

For standard coordinates, the symmetric faer KKT receives the physical
quadratic objective

```text
1/2 (eta_E / L^3) c^T K c
  + 1/2 sum_i precision_i (a_i^T x - target_i)^2
```

plus only the hard affine equalities and complete Cubic `Pi1` side conditions.
Soft scalar equalities therefore stay on the form-driven Symmetric KKT route;
they do not create QP slack variables or expose a backend control. Numerical
coordinate and Ruiz scaling remain derived reversible transforms and do not
define penalty, tolerance, or objective meaning.

## Recover and Verify

Recovery independently evaluates every hard and soft canonical expression in
physical units, reconstructs Cubic FieldEnergy with the configured `eta_E`,
computes each soft loss, and verifies

```text
total objective = 1/2 FieldEnergy + sum(soft losses).
```

The physical total is compared with a separately reconstructed standard-form
objective under the fixed `1e-11` round-trip envelope. Invalid recovery maps,
provenance association damage, and objective round-trip damage are distinct
`RecoveryVerificationReason` values. All reject the candidate before
`SolvedModel` construction.

## Evidence

`tests/public_soft_field_values.rs` proves checked scalar boundaries,
repairable build failure, public end-to-end fitting and querying, typed
standard-deviation/penalty equivalence, independent duplicate evidence, and
length/field-unit covariance of residuals and the complete objective.

`src/cubic_equality.rs` contains the narrow canonical recovery corruption test
`damaged_soft_objective_round_trip_is_rejected_without_a_model`, alongside the
existing damaged coordinate and provenance tests. The cumulative v0.1.0
public, canonical, failure and query suites remain unchanged and passing.
