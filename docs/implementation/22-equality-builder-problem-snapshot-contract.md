# Issue 22: Equality builder and problem snapshot input contract

Issue: [#22](https://github.com/qingsongyukuai/GeoRBF/issues/22)

Evidence seams: T03, T04

Requirements: PAPI-006–PAPI-009, PAPI-017; DOM-003, DOM-020–DOM-022;
IR-008

## Contract closure

This ticket closes the v0.1.0 public input contract assembled by issues #19–#21.
Checked geometry, frame, metric, observation, group, and gauge constructors keep
invalid local values outside `ProblemBuilder`. The sealed `ProblemInput` boundary
means every value reaching `ProblemBuilder::add` is already locally valid. An add
with a duplicate `SourceId` or `GroupId` rejects before mutating any identity set
or input collection, including a complete group containing several otherwise-new
members.

Forward references from additive gauges to a `GroupId` remain valid during the
mutable phase. `ProblemBuilder::build` resolves all such references together,
sorts missing-reference evidence by stable caller identity, and reports it with
independent problem-configuration errors. `BuildFailure::into_builder` returns
the original unsorted mutable input state so callers can add the missing groups,
repair configuration, and retry without reconstructing accepted inputs.

A successful build sorts observations by `SourceId` and semantic groups by
`GroupId` before creating the owning `Arc`-backed snapshot. The snapshot retains
the caller's frame, length and field units, global metric, input values, group
kind and membership, source/group identities, and fit configuration. It is
immutable, cheaply cloned, and `Send + Sync`; concurrent fits read the same
canonicalized snapshot without shared mutation.

Exact duplicate hard facts retain a separate source relation and public recovery
assessment for every caller `SourceId`, while canonical lowering may share their
one exact equality. The canonical equality count therefore reflects exact
deduplication without erasing provenance. No tolerance-coordinate merge,
epsilon perturbation, jitter, or soft-evidence merge is introduced.

## Property evidence

`tests/public_problem_contract.rs` adds a test-only `proptest` harness at the T03
and T04 public seams. It checks:

- arbitrary floating-point bit patterns across points, vectors, field values,
  gauges, tangent normalization, frames, and anisotropy metrics without panics;
- atomic duplicate `SourceId`, duplicate `GroupId`, duplicate group-member, and
  cross-group source rejection;
- empty group rejection, complete-group immutability, deterministic aggregation
  of multiple dangling references and resource errors, and repair after failure;
- randomized top-level input permutations with identical snapshot metadata,
  problem size, hard-relation provenance, canonical equality count, accepted
  objective, and model sample;
- preservation of both sources for an exact duplicate hard fact; and
- concurrent read-only fitting from cloned snapshots.

`proptest` is a dev dependency with only its `std` feature. The production crate
dependency and feature envelope is unchanged.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --test public_problem_contract
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```
