# Issue 16: production equality spine

Issue: [#16](https://github.com/qingsongyukuai/GeoRBF/issues/16)

Evidence seams: T14, T15

Requirements: PAPI-002, NUM-001, NUM-002, NUM-014, VAL-013, VAL-015

## Production admission verdict

`faer = 0.24.4` is admitted for the product-internal, sequential dense
symmetric-indefinite KKT route with only its `linalg` and `std` features. The
crate uses faer's low-level LBLT interface so every factor and solve allocation
is visible to GeoRBF's preflight plan. No other backend is admitted in v0.1.0.
GeoRBF resolves and owns the `georbf-v1` policy, `PartialDiag` pivoting,
64-column LBLT block, parallelism threshold, sequential execution, and `1e-11`
backend-standard-form backward-error limit. BackendFingerprint records those
resolved faer settings separately from the numerical-policy identity.

The previously ambiguous factor workspace is obtained exactly from the public
`cholesky_in_place_scratch::<usize, f64>` API with `Par::Seq` and faer's locked
default LBLT parameters. The plan separately records the returned allocation
size and alignment. It obtains solve scratch from the corresponding public
`solve_in_place_scratch` API and uses the larger of the sequential workspaces
in its conservative peak.

For the 3×3 manufactured KKT, factor and solve scratch are each 64 bytes with
64-byte alignment. The deterministic T14 boundary case accepts KKT dimension
23,151 and rejects adjacent dimension 23,152 at 8,589,948,160 planned bytes,
13,568 bytes over the 8 GiB limit. The factor workspace at the rejected shape
is 12,039,040 bytes. This test computes only layouts; it performs no large
allocation and invokes no backend. A separate `usize::MAX + usize::MAX` case
returns typed arithmetic-overflow evidence.

## Peak plan

Checked arithmetic covers canonical scalars; the dense Equality Hessian,
Jacobian, and right-hand sides; assembled KKT; factor matrix, subdiagonal, and
permutations; exact faer workspace; recovery buffers; and report storage. The
conservative peak holds all product components live and the larger of the
sequential factor/solve scratch layouts. A peak over 8 GiB or any failed checked
operation returns structured Capacity Exceeded evidence before assembly or a
backend call.

The Equality route creates neither CSC storage nor query scratch, so issue #16
does not add placeholder budgets for those future capabilities. It also makes
no reference-runner performance claim. Those parts of the broader NUM-014/T14
seams become applicable only when their end-to-end product routes are admitted.

## Minimal consumer evidence

The crate-internal consumer assembles and factors the manufactured system

```text
[2 0 1] [ 0.5]   [2]
[0 2 1] [-0.5] = [0]
[1 1 0] [ 1.0]   [0]
```

and checks the candidate against the independent analytic literal shown above.
It computes normalized backward error only as a dimensionless backend-standard-
form contract check; it does not claim physical-unit canonical acceptance or
produce a solved domain model. Its evidence owns the candidate, equality
multiplier, numerical-policy identity, verification scope, capacity plan,
factor/solve workspace, and a BackendFingerprint containing faer version,
features, resolved settings, scratch API, target, and requested/actual
sequential thread count.

## Replay

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
python3 scripts/audit.py
cargo tree --locked --target all -e features
cargo package --locked
```

The product workflow repeats the Rust 1.85.0 check, issue #16 behavior tests,
fail-closed dependency/license audit, feature graph, and packaging check on all
five native target families required by NUM-002. The audit pins the reviewed
lockfile identity, rejects any unreviewed build script, and then checks active
packages, features, native links, and licenses on each target.
