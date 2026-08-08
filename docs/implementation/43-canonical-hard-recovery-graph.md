# Issue 43: exact canonical hard-relation recovery graph

Issue: [#43](https://github.com/qingsongyukuai/GeoRBF/issues/43)

The `CanonicalCubicSolverForm` now owns one recovery graph for hard equalities.
Its columns are canonical atoms: physical functional dimension, exact binary64
support, value or gradient component, and semantic-latent identity. Consequently two
nearby supports, including observations with the same direction, remain
different functionals without a merge radius, snapping, tolerance-based
deduplication, or geometric clustering.

Compression-eligible solver rows are considered in canonical order. Elimination may propose
a dependency, but admission requires an exact reconstruction of every canonical
coefficient and the target. Verification accumulates binary64 products as
arbitrary-width signed integers, so an equality produced only by rounded
floating-point addition is rejected. Failure to prove the complete affine
identity retains the row as a solver constraint. Pre-existing relations whose
domain lowering explicitly marks them as physical-verification-only (for
example, an additional additive-gauge convention) remain verification-only;
the graph labels those paths as not having a complete affine reconstruction,
so they cannot be mistaken for issue-43 compression evidence.

The graph records the retained solver-row indices and recovery coefficients for
each canonical row. It also records a separate path, original target, and
relation-to-canonical sign for every original hard source provenance, including
its `SourceId`, semantic role, residual identity, and canonical-row association.
Exact duplicates gathered during lowering thus keep all caller-owned paths even
when opposite normalized signs share one canonical relation.

Equality KKT and Convex QP both iterate the same retained-row list. Each adapter
verifies the graph before accepting provenance, while physical recovery still
evaluates every canonical relation and every source report through the original
canonical problem. QP capacity planning uses the retained equality count. This
issue intentionally supersedes issue 38's temporary promise that adapter row
counts remain unchanged; activating the repository-wide `georbf-v2` policy
identity remains assigned to issue #49.

The symbolic rows and elimination basis are sparse in canonical atoms. This
keeps unrelated Dense Hermite supports out of one another's elimination work
and avoids materializing a canonical-row-by-global-atom dense matrix.

Regression coverage includes exact duplicates, a non-trivial consistent linear
dependency on both solver routes, a target that is only equal after binary64
rounding, distinct support atoms, equal gradients at close but different
supports through the public KKT and QP workflows, and actual directed normals
with the same direction at different supports.
