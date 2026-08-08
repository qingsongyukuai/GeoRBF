# Issue 42: Bounded polynomial and quotient precision rescue

Issue: [#42](https://github.com/qingsongyukuai/GeoRBF/issues/42)

Primary evidence seam: the public `ProblemBuilder` through snapshot, fit,
diagnostics, and model-query workflow. Cancellation and modal-sign conformance
remain behind the crate-private precision seam established by issues 40 and 41.

## Targeted canonical recomputation

An f64 polynomial decision in the rank gray band is recomputed from the
canonical functional terms as a four-mode double-double Gram problem. For a
quotient pivot whose outward-rounded interval crosses zero, deterministic
symmetric diagonal pivoting first factors every reliably positive mode. Only
the unresolved trailing Schur block is upgraded. Its basis directions are
expanded through the implicit Householder complement, and their entries are
recomputed from the original canonical Cubic generalized-functional pairings;
the rounded quotient Gram block is not promoted and decomposed again.

The pure-Rust arithmetic supplied by issue 41 records 106 precision bits. The
isolated block is never truncated: up to 64 modes are upgraded in full, while
65 or more produce Numerical Decision Gray Zone evidence immediately.
Canonical-pairing absolute operation scales seed explicit error radii, which
are propagated through Schur products, division, and square root. A sign is
accepted only when the resulting outward interval excludes zero.

## Reattachment and certificates

Strictly positive rescued modes are attached to the stable symmetric-pivoted
prefix. The basis retains the permutation as part of its private reversible
coordinate map. The resulting f64 basis must again pass the complete LLT
backward, positive-pivot, FieldEnergy identity, complete-Pi1 side-condition,
solver recovery, and canonical-response round-trip certificates before it can
reach solver-form assembly. No ridge, Gram jitter, or mode truncation is used.

Successful analysis exposes the upgraded start, complete mode count, precision,
and conclusion for polynomial and quotient rescue. Failure reports preserve the
same evidence. An upgraded exact zero is accepted as rank evidence only when
the canonical functional combination reconstructs to algebraic zero. Negative
curvature remains a non-positive representation failure. An interval that
still meets zero, and an over-capacity rescue, are Numerical Decision Gray Zone
rather than rank deficiency.

## Conformance coverage

The crate-private corpus checks the independent oracle's small-positive,
algebraic-zero, and negative Schur conclusions; a cancellation-scale positive
quotient mode is reattached and reruns all basis certificates. Separate tests
pin the 64/65 boundary without prefix truncation and verify that zero, negative,
and unresolved conclusions remain distinct.
