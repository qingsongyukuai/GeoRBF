# Issue 39: Implicit Householder complete Pi1 quotient

Issue: [#39](https://github.com/qingsongyukuai/GeoRBF/issues/39)

Primary evidence seam: the existing public `ProblemBuilder` to snapshot to
`fit` to diagnostics/model-query workflow. The representation diagnostics also
expose construction counts and numerical defects without exposing matrices or
backend storage.

## Complete polynomial pairing

`CubicRepresentation::build` assembles the pairing of every retained fitting
functional with the complete `Pi1 = span{1, x, y, z}` basis. Exact functional
identity remains the established representer-span identity rule; no coordinate
radius, support spacing, or direction/normal heuristic merges distinct fitting
functionals. The reported representer, polynomial, and quotient dimensions are
taken from the corresponding constructions. A full-rank build therefore
reports `quotient_dimension = fitting_functional_count - 4`.

The unpivoted Householder QR of the `n by 4` polynomial pairing represents both
the rank-four polynomial image and its orthogonal complement. Only the four
polynomial-image columns are realized for the orthogonality certificate. The
quotient image `Q2` remains implicit in response projection and recovery.

## Direct trailing-block congruence

The quotient Gram matrix is formed by applying the stored Householder sequence
to the complete kernel pairing from the left and right, then copying the
trailing quotient block of `Q^T K Q`. This takes two full-matrix Householder
applications with four reflectors, so the congruence work is `O(n^2 p)` for
`p = 4`.

The construction reports the actual reflector count and increments the
congruence-pass count around the two full-matrix applications. Those counts
give a directly checkable `2p` operation shape. The implementation contains no
quotient-column construction loop or per-column kernel matrix-vector path; the
previous dense `Q2` columns and their full kernel multiplications have been
removed.

## Numerical evidence

Every successful construction records:

- the orthogonality defect of the realized rank-four Householder polynomial
  image;
- the trailing defect of `Q^T P`, which verifies the complete `Pi1` side
  condition without realizing `Q2`;
- a quotient-coordinate round trip for a canonical fitting-functional response;
- the actual representer, polynomial, and quotient dimensions; and
- the construction counts described above.

The numerical policy owns the construction limits: representation construction
rejects orthogonality and canonical-response round-trip defects above `1e-11`,
while the existing side-condition gate is stricter at `1e-12`. The public
regression independently requires all three defects to be at most `1e-11`.
Existing Equality KKT and Convex QP public workflows continue to consume the
same solver-independent canonical form, with the QP path now consuming the
directly transformed trailing block.
