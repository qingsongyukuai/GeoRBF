# Issue 38: Cubic Representation and Canonical Solver Form seam

Issue: [#38](https://github.com/qingsongyukuai/GeoRBF/issues/38)

Primary evidence seam: the existing public `ProblemBuilder` → snapshot →
`fit` → diagnostics/model-query workflow.

## Deep representation Interface

The crate-private `CubicRepresentation` is the Cubic Quotient Representation
Module. Its caller-facing responsibilities are deliberately limited to:

- `build`: construct the current complete Cubic representation from all
  fitting functionals and the resolved metric/FieldEnergy normalization, and
  initialize the solver-independent field portion of the canonical form;
- `response`: map one canonical field functional into the representation's
  retained solver-coordinate realizations and complete `Pi1` response;
- `recover`: map solver field coordinates back to canonical representer
  coefficients and the recovered physical query field, including coordinate,
  side-condition, coefficient, polynomial, and FieldEnergy round trips.

Householder storage, kernel/polynomial pairings, reduced energy coordinates,
faer types, and physical coefficient conversion remain implementation details.
No public matrix, coefficient, basis, backend, or adapter interface is added.

## One canonical solver form

`CanonicalCubicSolverForm` is the sole solver-independent assembly of field
responses, semantic-latent columns, FieldEnergy, hard equalities, affine
inequalities, soft objectives, and source recovery metadata. It records every
canonical row's index, provenance, residual/derived identities, target,
dimension, and participation before either solver route is selected.

The form retains the two numerical coordinate realizations already used by
v0.2.0 so this boundary-only change does not alter public behavior:

- Equality KKT consumes the existing standard representer coordinates and
  complete four-row `Pi1` side condition.
- Convex QP consumes the existing Householder quotient coordinates.

Both adapters now only realize their scalar dense backend inputs from the same
form. They no longer compute functional responses or independently lower
hard/soft semantics. Provenance verification compares each adapter realization
back to the immutable canonical form rather than reassembling a second response.
The form obtains every canonical functional's field coordinates only through
`response`; neither it nor either adapter can access representation matrices or
basis objects. The field-energy and side-condition coefficients initialized by
`build` are owned immediately by the canonical form rather than exposed through
a representation matrix getter or a non-functional response variant.
This seam lets issues #39 and #40 replace the representation coordinates without
changing canonical hard/soft ownership or either backend adapter.

## Compatibility and verification

The existing field dimensions, KKT/QP sizes, attempt policies, successful
models, structured failures, canonical physical acceptance envelope, soft-loss
associations, diagnostics, and query path are unchanged. Representation
recovery is used by both Equality and QP, while all-source physical verification
continues to evaluate the recovered field through the canonical problem.

Replay under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```
