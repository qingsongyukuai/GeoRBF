# Issue 24: immutable solved-model logical batch queries

Issue: [#24](https://github.com/qingsongyukuai/GeoRBF/issues/24)

Evidence seam: T09

Requirements: PAPI-014, PAPI-015, PAPI-018; IR-013; NUM-014; VAL-006

## Owning immutable model contract

`SolvedModel` owns one `Arc` containing the recovered Cubic field, the complete
owning `ProblemSnapshot`, and recovered shared-level values. Cloning a model
therefore clones the `Arc`, while every query is read-only and the contained
state is `Send + Sync`. The retained snapshot exposes the resolved Cubic Kernel
Contract, `NumericalPolicyId`, input coordinate frame, length unit, and field
unit. Recovered horizon and shared-level values remain addressable by stable
`GroupId`. Model debug output and public accessors do not expose representer or
polynomial coefficients, center ordering, matrices, or backend variables.

## One single/batch evaluation meaning

`SolvedModel::evaluate_batch(&[Point3])` returns an owning `Vec<FieldSample>` in
the caller's exact input order. Empty batches return an empty vector and
repeated points remain repeated. Every batch element calls the same
`SolvedModel::evaluate` path as a single-point query, which in turn evaluates
the recovered field through the Cubic generalized-functional pairing and
returns one coherent field value and complete three-component gradient in the
input frame. The batch path does not add query functionals to the fitted
representer span and cannot mutate the model or its snapshot.

The result vector is local until all points succeed. A non-finite recovered
observable therefore drops the entire staged result and returns `QueryError`
with `NonFiniteResult` and the first failing logical-batch index. A single-point
failure has no batch index. Finite large coordinates are accepted whenever all
recovered observables remain finite; overflow is a typed query failure rather
than a partial result or panic.

## Checked scratch and streaming centers

Atomic result staging and the temporary internal-sample chunk are planned
before allocation with checked arithmetic and a combined hard 256 MiB
query-scratch limit. A representable over-limit plan reports both planned bytes
and the limit; arithmetic overflow reports the same typed `CapacityExceeded`
reason with no fabricated planned size. Allocation and point evaluation begin
only after the plan succeeds.

The private query plan bounds the temporary chunk of recovered value-gradient
samples, but chunk size is not public API and does not alter order or arithmetic.
Within each chunk, every point streams across the recovered representers and
accumulates only its value and three gradient components; no
logical-batch-by-center value, gradient, or jet matrix is materialized.
Consequently center count affects work, not batch scratch. The plan-level
boundary test admits 100,000 returned value-gradient samples, checks the exact
256 MiB boundary including a one-sample internal chunk, and rejects the first
over-limit plan before allocation.

## Validation evidence

Public T09 coverage verifies exact single/batch equivalence and input order;
empty, repeated, finite large, and overflow-producing coordinates; atomic first
index failure; multiple batch shapes spanning internal chunk boundaries; and a
100,000-location smoke with a reduced center set. Concurrent threads mix single
and batch calls through several cheap model clones and compare only canonical
field observables. Automatic and exact-one fit resource plans retain the same
model contract and agree within VAL-006's
`1e-12 * field_scale + 1e-11 * sample_reference_scale` envelope.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --test public_model_queries
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```
