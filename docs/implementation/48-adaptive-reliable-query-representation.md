# Issue 48: adaptive reliable query representation

Issue: [#48](https://github.com/qingsongyukuai/GeoRBF/issues/48)

Primary evidence seam: the public `ProblemBuilder` → `fit` →
`RepresentationEvidence` → single/batch model-query workflow.

## Verified query recovery

A successful recovery now retains one immutable query expansion containing the
physical generalized-RBF coefficients and the centered/scaled complete Pi1
coordinates from the stable solver representation. The expansion becomes
queryable only after all-source recovery, the physical Pi1 side condition,
polynomial and coefficient basis round trips, and the FieldEnergy round trip
have passed. `RepresentationEvidence::verified_query_representation` records
those conclusions and the maximum basis round-trip error.

Keeping the centered Pi1 coordinates prevents a large physical constant from
destroying low-order field information at translated query locations. Recovery
still retains physical polynomial coefficients for canonical verification; the
query representation is a stable evaluation of that same recovered field, not
a second fit or a model choice.

## Adaptive deterministic accumulation

Value and all three physical gradient components use the same ordered Neumaier
accumulation over polynomial and generalized-RBF contributions. Each component
is accepted only when its rounding bound fits
`1e-12 * field_scale + 1e-11 * abs(sample_component)`.

When that evidence is insufficient, the query streams the same polynomial,
coefficients, representers, and Cubic jets once more with the existing
pure-Rust double-double arithmetic. No coefficient changes, quotient solve,
matrix materialization, fit, or mutable state occurs. If the bounded retry
still cannot certify every component, the query returns the non-exhaustive
`QueryErrorReason::NumericalIndeterminate` variant. Genuine non-cancelling
overflow remains `NonFiniteResult`.

Single queries and every batch chunk call this exact evaluation path. Batch
results remain staged until all points succeed, retain input order, and attach
the first failing logical index to either numerical failure. Work is O(n) per
point and scratch remains independent of representer count.

## Regression evidence

The public model-query suite checks the successful representation-evidence
bundle, a five-support Cubic cancellation field against independently computed
high-precision value and gradient literals, structured indeterminate single and
atomic-batch failures, and the existing chunking, ordering, capacity,
concurrency, and non-finite behavior.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --test public_model_queries
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
```
