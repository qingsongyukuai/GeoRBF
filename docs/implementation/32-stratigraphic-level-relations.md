# Issue 32: Stratigraphic Age and Field Level Relations

Issue: [#32](https://github.com/qingsongyukuai/GeoRBF/issues/32)

Primary evidence seams: PAPI-008–PAPI-010; DOM-002–DOM-003, DOM-012–DOM-013,
DOM-015, DOM-017, DOM-019–DOM-020; IR-004–IR-005, IR-008–IR-010, IR-014;
NUM-009–NUM-010; DIA-002, DIA-006–DIA-009; VAL-004, VAL-008, VAL-015

## Checked public semantics

`StratigraphicFieldDirection` is an explicit, single-assignment problem
configuration with `TowardYounger` and `TowardOlder` values. A snapshot that
contains `YoungerThan` or `OlderThan` cannot be built without it. Repeating the
same direction or attempting to replace it is atomically rejected; neither
age nor field-value orientation is inferred from coordinates, handedness,
axis labels, input order, names, or `GroupId` ordering.

`MinimumFieldSeparation` accepts only finite, strictly positive values in the
problem's field-value units. `YoungerThan` and `OlderThan` always carry one.
`FieldLevelOrder` is a separate non-strict relation whose declared lower and
upper shared levels may recover equal values. Each relation offers only hard,
quadratic-violation, and linear-violation constructors; soft relations require
the snapshot's explicit `FieldEnergyNormalization`.

All three inputs retain a caller-owned `SourceId` and two stable `GroupId`
roles. They can be inserted before their groups. Build aggregates unresolved
references in deterministic `SourceId`/`GroupId` order and retains the builder
for repair. A referenced singleton shared level is informative because its
semantic latent now participates in a relation.

## Relation graph preflight

Hard relations normalize to directed difference edges in field-increasing
orientation. Age edges carry their positive minimum separation; Field Level
Order edges carry zero. Preflight rejects hard self-relations and every cycle
containing a strict edge, including direct reverse contradictions and chains
closed by a non-strict order. Stable evidence retains every constituent
`SourceId`, `GroupId`, semantic role, and the fact that no backend ran.

The same graph accumulates minimum separations along hard chains. When its end
groups have absolute representatives supplied through a gauge, a field-value
observation, or the existing shared-level equality forest, an incompatible
finite difference is proved before lowering reaches a backend. Soft edges do
not participate in conflict proofs and remain independent violation evidence.

## Canonical lowering and Cubic QP recovery

Every legal relation lowers into the one physical `CubicCanonicalProblem` as
an affine inequality over the two stable semantic latents. The canonical row
is `upper_field_level - lower_field_level >= required_difference`, where the
required difference is the checked minimum for age relations and exactly zero
for Field Level Order. Its coefficients sum to zero, so the relation remains
additive-gauge invariant. Any affine inequality selects the capability-driven
Clarabel QP route; hard duplicates share canonical mathematics while retaining
all source and both-group provenance, and soft duplicates retain independent
violation variables and objective terms.

Recover and Verify independently restores both shared values, their oriented
field separation, physical slack or violation, tolerance, active state, loss,
and complete relation provenance. `FitReport::shared_level_relations` is sorted
by `SourceId` and reports the original relation kind and group roles together
with the resolved field direction. `SolvedModel::shared_level_value` and normal
single or batch queries use the same accepted field and snapshot.

## Conformance evidence

`tests/public_stratigraphic_relations.rs` covers checked inputs, atomic
configuration, forward and dangling references, self/reverse/strict/accumulated
graph conflicts, both field directions, non-strict equality, hard and soft
duplicates, quadratic and linear violations, absolute gauge/value composition,
input and ID permutation, frame metadata independence, stable reporting,
queries, and certificate-validated general infeasibility.

Narrow contract tests inject QP scaling-recovery corruption and checked
capacity failure through a public relation snapshot. They require structured
failure evidence and no model. The issue 30/31 QP, recovery, provenance, and
cumulative public suites remain part of this capability's release boundary.
