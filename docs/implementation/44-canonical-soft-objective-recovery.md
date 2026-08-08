# Issue 44: canonical soft-objective recovery

Issue: [#44](https://github.com/qingsongyukuai/GeoRBF/issues/44)

The `CanonicalCubicSolverForm` now owns a soft recovery graph alongside its
hard recovery graph. Every canonical soft scalar residual is retained as a
solver row and is associated with its original provenance, residual identity,
objective block, and component position. The graph also snapshots the complete
canonical solver rows, including functional responses, targets, and physical
dimensions, plus the complete objective definitions, including loss, precision,
whitening, covariance-group identity, and block kind.

Soft rows deliberately remain one-for-one. Equal row spaces or identical
functionals do not prove objective identity because multiplicity, targets,
weights, and block-level losses remain independently meaningful. A future soft
compression therefore cannot become active merely by sharing the hard-row
elimination machinery: it must replace this identity map with an exact
all-candidate-field objective proof and a complete per-source recovery map.

Both Equality KKT and Convex QP reject their assembled forms unless the soft
graph still verifies. Verification requires canonical row order and complete
coverage; exact source provenance and residual identities; exact objective and
component associations; unchanged functional responses, targets, and physical
dimensions; and unchanged loss, precision, whitening, covariance, and block
semantics. Hard recovery remains a separate graph, so a hard relation
with the same functional cannot absorb or replace a soft objective row.

Physical recovery continues to evaluate every original soft equality in its
original units. Independent observations expose their own residual, configured
loss, and objective contribution. Covariance groups expose every member's
source-associated residual while retaining one canonical group-level objective;
they do not invent a non-unique per-member loss split.

Regression coverage checks graph rejection after objective or provenance
mutation, duplicate soft field-value observations with different weights, and
hard/soft coexistence on the same functional through both Equality KKT and
Convex QP. Existing covariance-group tests continue to verify group-level loss
and per-member residual recovery.
