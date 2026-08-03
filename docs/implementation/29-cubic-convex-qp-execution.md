# Issue 29: Cubic Convex QP execution seam

Issue: [#29](https://github.com/qingsongyukuai/GeoRBF/issues/29)

Primary evidence seams: T03–T07, T13, T15

Requirements: IR-001, IR-004–IR-012, IR-014; KER-006–KER-008; NUM-001,
NUM-003–NUM-015; DIA-001, DIA-003, DIA-005; VAL-005, VAL-013, VAL-015

## Capability-driven execution

The existing physical `CubicCanonicalProblem` is the only domain model. A
crate-private algebraic planner chooses the symmetric Equality KKT when its
relations are affine equalities and chooses Convex QP only when the same
canonical problem contains an affine inequality. The plan records form family
and relation counts but owns no faer, Clarabel, sparse-layout, or raw backend
type. This issue deliberately adds no public relation, solver selection,
matrix, coefficient, or backend API.

Both routes return the shared recovered field, semantic latent, provenance,
hard/soft relation, FieldEnergy, objective, side-condition, and canonical
acceptance concepts. The manufactured contract consumer adds a redundant upper
bound to an otherwise identical canonical problem. Its KKT and QP executions
agree on field coefficients, polynomial, semantic latent, FieldEnergy, soft and
total objectives, and value/gradient queries within the fixed validation
envelopes.

## Cubic QP reduction and realization

The QP route reuses the complete Cubic representation analysis. It computes the
full `Pi1`, classifies rank with the shared RRQR/SVD gray bands, constructs the
implicit Householder null-space map, and verifies null defect, affine
reproduction, reduced symmetry, and positive definiteness before Clarabel may
be entered. Negative curvature, rank loss, and a numerical gray zone remain
structured representation evidence; no ridge, jitter, mode deletion, kernel
substitution, or other hidden repair is permitted.

The solver-independent dense form contains only scalar vectors, canonical row
provenance, and semantic block metadata. It assembles

```text
1/2 z^T Q z + q^T z
subject to Aeq z = beq and Aineq z <= bineq,
```

where the field block is the normalized reduced Cubic energy and soft blocks
contribute `A^T P A` and `-A^T P target`. Only
`src/clarabel_backend.rs` converts this form to Clarabel CSC matrices, cones,
settings, and raw primal/dual/slack vectors. Clarabel is pinned to 0.11.1 with
default features disabled and only `serde` selected; qdldl is fixed as the
direct linear solver and one thread is required and checked.

## Capacity, scaling, and attempts

Checked capacity planning completes before representation/form allocation or
backend entry. Its conservative peak includes physical canonical storage,
dense QP assembly, upper-Hessian and constraint CSC realization, qdldl factor
workspace, recovery, and report storage. Arithmetic overflow and a peak above
8 GiB return structured evidence with both allocation and backend-entry flags
false.

GeoRBF applies exactly eight deterministic Ruiz-style rounds. Every row and
variable factor is an exact power of two, each round is clamped to exponent
`[-8, 8]`, cumulative exponents are clamped to `[-32, 32]`, and forward/inverse
maps are retained and verified. Clarabel's own equilibration is separately
fingerprinted and its factors are used only for a backend-internal round-trip
check; they never redefine physical tolerances or canonical semantics.

The immutable attempt plan is Standard then Robust, at most once each. Both
profiles fix the full settings fingerprint, qdldl, and one thread. A later
attempt changes only numerical settings: canonical tolerances, objective,
constraints, and QP form family remain unchanged. Every attempt records status,
settings, threads, iterations, backend reports, internal scaling, and
independently recomputed primal, dual, stationarity, complementarity, and gap
residuals. A backend status alone never accepts a candidate.

## Recover and Verify

Recovery reverses GeoRBF scaling and the Householder reduction, reconstructs
the physical field, polynomial, and semantic latents, and then independently
checks:

- backend standard-form rows and nonnegative inequality slack;
- physical hard-equality and affine-inequality violations by dimension;
- canonical row/provenance association;
- scaling, reduction, polynomial, coefficient, FieldEnergy, whitening, and
  objective round trips;
- Cubic side conditions, finite recovered quantities, and total objective
  `1/2 FieldEnergy + 1/2 sum(r^T P r)`.

Injected provenance, GeoRBF scaling, Householder recovery, backend-residual, and
objective corruptions each produce structured failure evidence and no model.
An injected unverified Standard attempt proves the single deterministic Robust
retry. A flattened manufactured `Pi1` proves a missing affine mode is reported
before Clarabel entry.

## Dependency and cumulative evidence

`scripts/audit.py` pins the exact LF-byte lockfile identity, production
dependencies, Clarabel/faer feature sets and versions, selected build scripts,
licenses, native links, and forbidden native/BLAS/LAPACK/SDP packages and
features. The production closure remains pure Rust and compatible with Rust
1.85. The cumulative tests, strict Clippy, rustdoc, packaging, and release build
remain the acceptance boundary; the public v0.1 Equality API and its behavior
are unchanged.
