# Issue 47: representation evidence and diagnosis semantics

Issue: [#47](https://github.com/qingsongyukuai/GeoRBF/issues/47)

Every `FitReport` now publishes one `RepresentationEvidence` bundle through
the same public fit-report seam on success and failure. The bundle identifies
the executed `georbf-v2` policy and records the failing and last completed
construction stages. Canonical hard and soft dimensions, participating
`SourceId`s, representers, complete Pi1 dimension and rank, quotient dimension,
retained/truncated modes, and recovered-source coverage retain their owning
construction's counts instead of borrowing a solver row count.

Completed construction exposes the existing implicit-Householder certificate
and the verified unregularized LLT bundle. Together these report Householder
orthogonality, canonical-response round trips, pivot intervals, normalized
backward error, FieldEnergy identity, Pi1 side condition, solver-coordinate
recovery, retained modes, and forbidden ridge/jitter/truncation flags. Bounded
polynomial or quotient rescue evidence is lifted into the common bundle. A
pre-backend failure retains its structured analysis failure there as well, so
the actual quotient dimension, pivot interval or rescue conclusion and all
sources that reached representation participation remain auditable.

`AllSourceRecoveryEvidence` is embedded after recovery. It independently lists
participating and recovered sources and proves whether the unregularized
canonical problem passed complete recovery. Problem regularization remains
unconditionally false. Backend factorization regularization is recorded per
attempt and stays separate: Equality KKT attempts disable it, while Clarabel
attempts transparently report whether their internal static or dynamic
stabilization was enabled. Such a setting does not imply canonical acceptance;
only Recover and Verify can set the unregularized recovery result.

Problem diagnosis no longer folds a proved representation negative-curvature
direction into generic numerical failure. It remains distinct from algebraic
rank deficiency, bounded-rescue gray zone, backend contract violation,
infeasibility, and recovery verification failure. Condition estimates,
spectrum summaries, KKT rank/inertia, and refinement remain attempt and risk
evidence; the representation diagnosis mapping consumes only canonical mode,
precision-rescue, curvature, and recovery evidence.
