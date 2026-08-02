# Issue 17: product-internal Cubic Equality core

Issue: [#17](https://github.com/qingsongyukuai/GeoRBF/issues/17)

Evidence seams: T05, T06, T13

Requirements: IR-002, IR-003, IR-010–IR-012; KER-001, KER-002,
KER-004–KER-009; NUM-006, NUM-010; VAL-001, VAL-003–VAL-005, VAL-013,
VAL-015

## Product boundary

The `georbf` crate now owns an end-to-end, crate-internal Cubic Equality path.
It accepts finite value/first-derivative scalar functionals, constructs the
Cubic generalized-functional representation, verifies the complete affine
polynomial space, assembles the symmetric Equality KKT through the production
`faer` spine admitted by issue #16, and recovers a queryable field. No public
functional, matrix, coefficient, polynomial, kernel-extension, or backend API
is added by this milestone.

The canonical functional normal form sorts exact supports, merges only exactly
equal supports, deletes zero coefficients and terms, canonicalizes signed zero,
and returns a typed `ZeroFunctional` failure when cancellation removes the
whole functional. Output dimension is stored on the functional. Typed source,
optional group, relation, and typed semantic-role-path provenance are stored only on each
usage edge, so structurally equal functionals can be interned without losing
relation provenance.

## Cubic contract and oracle adoption

For `d=x-y`, `v=Md`, and `r=sqrt(d^T M d)`, the implementation owns

```text
k(x,y)       = r^3
grad_x k     = 3 r v
grad_y k     = -3 r v
grad_xgrad_y = -3 (r M + vv^T/r)
```

The exact coincident-support branch returns canonical positive zero for the
value, both first jets, and the complete mixed jet. Nonzero displacement is
scaled before evaluating the quadratic form, so an underflowed squared distance
cannot be mistaken for the origin. `GlobalAnisotropyMetric` accepts only finite,
exactly symmetric, positive-definite 3×3 metrics with a finite determinant of
one under the owned policy; anisotropy enters only through `M` in the formulas
above.

The T05 contract tests consume byte-identical adoptions of issue #15's
independently generated 120-digit fixtures under
`validation/oracle/cubic-v1`. Tests verify the source declarations, fixture
bytes, stable CaseIds, provenance, output hashes, precision/rounding metadata,
and exact hexadecimal `f64` values before comparing the product value, first
jets, mixed jet, origin branch, generalized-functional contraction, and affine
observations. Analytic regressions additionally cover exchange/derivative signs,
frame/metric covariance, a four-term difference, and normal/tangent
functionals. The independent Python generator and verifier remain disposable
outside the product crate; only their generated declarations, fixtures, and
source manifest are adopted.

## Complete Pi1 and Equality assembly

The representation span is formed only from functionals used by the hard
fitting equalities. A deterministic coordinate center and positive common
length define the numerical polynomial coordinates. The plan retains the full
`Pi1 = span{1,x,y,z}` pairing `P` and an explicit reversible map between these
coordinates and physical polynomial coefficients.

Before the Equality solve, the core uses SVD with the owned `georbf-v1`
accept/reject bands to require `rank(P)=4`. A sequential Householder QR
materializes only the small T06 evidence pairing `T^T K T`, verifies
`range(T)=ker(P^T)`, checks its normalized symmetry defect against `256 u d`,
and requires strict positive definiteness. Four augmented interpolation solves
then verify affine reproduction. A rank-three manufactured case is rejected
before the field solve and records that neither a solver nor hidden
regularization was used.

For the ten-functional manufactured case, the physical primal contains ten
field coefficients followed by four explicit normalized polynomial
coefficients. The equality Jacobian contains four `P^T c=0` side-condition rows
followed by all ten hard generalized-functional rows `[K P]`. The Hessian keeps
the complete `K` block and a zero polynomial block. Issue #16 therefore solves
the resulting 28×28 symmetric augmented KKT; there is no ridge, jitter,
automatic anchor, polynomial deletion, or kernel fallback in the construction.

## Recovery evidence

The T06 manufactured targets are independent 120-digit Decimal evaluations of
the issue #14 declaration. Recovery checks the physical field coefficients,
the reversible physical polynomial, `P^T c`, every hard equality, the Cubic
quotient FieldEnergy, and canonical value/gradient observables. Hard-equality
recovery evaluates the recovered physical field through each canonical
functional instead of reusing the assembled K/P rows, and returns its typed
dimension and complete usage provenance with the value, target, and residual.
Field-value and field-value-per-length residual maxima remain separate and are
checked against independently owned physical-unit acceptance fields; they are
never folded into a dimensionally invalid raw maximum.
The side condition is reconstructed from the canonical functionals, the Cubic
energy pairing is recomputed from the recovered representers, and the physical
polynomial is mapped back to normalized coordinates for an explicit `1e-11`
round-trip check.
The all-hard FieldEnergy normalization resolves to one through a dedicated
Cubic unit-covariant type; its `eta_E' = s^3 eta_E / t^2` transformation is
tested against the corresponding native-energy transformation. Acceptance uses
the frozen `1e-11` backend backward-error limit, `1e-10` side-condition limit,
and `1e-8` canonical recovery envelope from the owned numerical policy.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```

The product workflow runs the accumulated issue #16 and #17 behavior suite on
all five admitted native targets. The lockfile and production dependency graph
are unchanged: `faer = 0.24.4` remains the only product dependency and uses the
same pure-Rust `linalg,std` feature envelope.
