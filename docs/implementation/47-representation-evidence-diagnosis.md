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
Checks inside the factorization certificate are optional until executed, so an
early LLT rejection cannot masquerade as an exact zero-error energy, side,
recovery, or response check. Response-assembly rejection retains the completed
Householder and LLT certificates and identifies response assembly as its stage.

`AllSourceRecoveryEvidence` is embedded after recovery. It independently lists
participating and recovered sources and proves whether the unregularized
canonical problem passed complete recovery. The representation construction
directly reports whether problem regularization was applied; the current policy
reports false. Backend factorization regularization is recorded per attempt and
stays separate. Equality KKT attempts report both mechanisms disabled. Clarabel
reports static regularization as applied when configured and explicitly marks
dynamic regularization as enabled-but-not-reported when its adapter cannot
observe whether that mechanism fired. Such a backend status does not imply
canonical acceptance; only Recover and Verify can set the unregularized
recovery result.

Problem diagnosis no longer folds a proved representation negative-curvature
direction into generic numerical failure. It remains distinct from algebraic
rank deficiency, bounded-rescue gray zone, backend contract violation,
infeasibility, and recovery verification failure. Condition estimates,
spectrum summaries, KKT rank/inertia, and refinement remain attempt and risk
evidence; the representation diagnosis mapping consumes only canonical mode,
precision-rescue, curvature, and recovery evidence.
