# Issue 46: all-source recover and verify

Issue: [#46](https://github.com/qingsongyukuai/GeoRBF/issues/46)

Every candidate-producing backend now crosses one solver-independent all-source
verification boundary after recovery to the canonical physical field. Equality
KKT and Convex QP pass their recovered hard relations, soft residuals, objective
blocks, and route-specific affine relation count to the same verifier owned by
the canonical solver form.

The verifier expands every hard recovery edge back to its original relation
orientation and applies that relation's Canonical Physical Acceptance Envelope.
Consequently an exact duplicate or compressed dependency cannot inherit success
merely because its retained solver row passed. Soft verification preserves the
complete canonical component indices, configured loss, covariance association,
and residual-block kind, so multiplicity and objective ownership remain intact.

Accepted fit reports expose `AllSourceRecoveryEvidence`. Its canonical hard and
soft relation counts, sorted participating and recovered `SourceId` sets,
representer count, solver-independent relation-row count, and recovery-edge
count come from their owning constructions. None is inferred from another
dimension. Acceptance requires exact equality between independently collected
participation and recovery source sets, verified hard and soft recovery graphs,
complete recovered relation counts, per-source hard acceptance, and unchanged
soft objective associations.

A missing source edge, damaged recovery graph, rejected source-level hard
relation, or changed objective association adds a recovery-verification reason
and prevents model construction. Existing round-trip, finiteness, side
condition, FieldEnergy, whitening, and objective-identity gates remain part of
the complete candidate boundary.

Regression coverage exercises consistent hard duplicates, duplicate soft
observations with distinct weights, hard and soft observations sharing one
functional, and close supports through both public KKT and QP fit/report paths.
