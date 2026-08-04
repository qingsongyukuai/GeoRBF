# Issue 36: Convex Relations preview release

Issue: [#36](https://github.com/qingsongyukuai/GeoRBF/issues/36)

Evidence seams: T01–T17 within the v0.2.0 Cubic KKT/QP scope

Requirements: PAPI-001–PAPI-015, PAPI-017–PAPI-019; DOM-001–DOM-007,
DOM-009–DOM-022; IR-001–IR-014; KER-001, KER-004–KER-009;
NUM-001–NUM-015; DIA-001–DIA-009; VAL-001–VAL-015.

## Public Convex Relations preview

The runnable `examples/convex_relations.rs` is the cumulative public acceptance
seam. It constructs hard and soft Field Value, complete Gradient, and Tangent
observations; a named Covariance Group; hard and soft Field Value Bounds and
Directional Derivative Intervals; ordered Horizons with a minimum separation;
Field Level Order; Field Separation Intervals; Point-to-Level-Set relations;
and both Directed and explicitly resolved Axial Normal observations.

All inputs enter one immutable `ProblemSnapshot`. Algebraic capability alone
selects the shared Cubic `Pi1` KKT or QP form. A successful candidate must pass
backend-contract checks, reversible scaling/reduction, physical canonical
recovery, objective reconstruction, provenance verification, and the unchanged
acceptance envelope before the example receives a `SolvedModel`. It then audits
typed relation assessments and shared values and compares ordered batch queries
with the corresponding single value-gradient queries.

The cumulative public corpus retains every v0.1.0 case and the mandatory,
pairwise, and higher-order cases delivered by issues 27–35. In particular it
covers covariance/vector/frame interactions; Cubic/`Pi1`/gauge combinations;
interval/order/shared-latent conflicts; normal direction and slope; retry,
certificate, recovery, and provenance behavior; and failure outcomes that
never contain a partial model.

## Capacity and smoke evidence

Checked planning continues to enforce an 8 GiB fit peak, a 256 MiB query
scratch limit, and the 100,000-location logical batch contract. The v0.2.0
release seam additionally runs a 512-scalar-constraint Cubic QP with 10,000
ordered queries. Five hard values establish the manufactured field and 507
distinct finite-support bounds keep all 512 scalar constraints in the
pre-presolve Convex QP quantity envelope.

The smoke is a regression tripwire, not a wall-clock SLA. The manifest-pinned
bare-metal runner, frozen absolute workloads, nine-process statistics, and a
formal `ValidationProfileId` remain later release-evidence work. This preview
therefore honors VAL-009–VAL-012 by keeping ordinary CI timings non-normative
and by making no formal performance or workload-stability claim.

## Traceability audit

`validation/v0.2.0/traceability.json` contains exactly one record for each of
the 99 scoped requirements. Every record names one unique behavior, one unique
public API or Canonical IR path, a concrete test/oracle evidence set, and an
implementation or release document. `scripts/release_traceability.py` rejects
missing, duplicated, unexpected, dangling, or unused records and reference
sets. The repository audit also verifies the independent oracle mirror byte for
byte, scans product source for unsupported placeholders, checks package
metadata and contents, and requires the public example and release workflow.

The v0.1.0 trace remains tracked as historical evidence. It is not used to
weaken or replace the cumulative v0.2.0 map.

## Compatibility and exclusions

This is a `0.x` Rust integration preview. Patch releases will not silently
change the published mathematical meaning, diagnosis categories, stable
identity semantics, or `georbf-v1` numerical policy. Source compatibility may
change in a later `0.x` minor release and the final v1 API is not frozen.

The release contains no Directed Normal Cone, Angle, Cone Violation, vector
linear-norm epigraph, SOC/SOCP product route, Gaussian, Inverse Multiquadric,
Wendland C2, advanced functional, reference-runner SLA, or v0.3.0+ placeholder.
It does not publish to crates.io.

## Release procedure

1. Run local typechecking and targeted tests throughout implementation, then
   the full test suite exactly once at the end.
2. Commit the reviewed release artifacts and push the accepted commit.
3. Dispatch the `release` profile of `.github/workflows/product-v0.2.yml` for
   that exact commit and require every job to pass with zero waiver.
4. Create an annotated `v0.2.0` tag targeting the accepted commit and push it.
5. Require the tag-triggered release profile to pass with zero waiver.
6. Create the GitHub release from `RELEASE_NOTES.md`, verify that its target and
   the annotated tag resolve to the accepted commit, then close issue #36 with
   links to both accepted workflow runs.
