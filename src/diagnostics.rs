//! Typed public fit diagnoses and backend-attempt evidence.

pub use crate::cubic_equality::RecoveryVerificationFailureReason as RecoveryVerificationReason;
pub use crate::functional::SemanticRolePath;
use crate::functional::{GroupId, SourceId};
use crate::geometry::Vector3;
pub use crate::numerical::NumericalPolicyId;

/// GeoRBF's semantic conclusion for a failed fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProblemDiagnosis {
    /// One or more caller inputs require an explicit semantic resolution.
    UnresolvedSemantics,
    /// The supplied problem could not form a valid fitting problem.
    InvalidProblem,
    /// Two locally comparable hard inputs have incompatible exact targets.
    DirectInputConflict,
    /// Shift-invariant relations did not select an absolute field representative.
    UnidentifiedAdditiveGauge,
    /// A one-member shared level set did not constrain or connect its latent.
    UninformativeSharedLevelSet,
    /// The observations did not identify every observable field mode.
    UnidentifiedFieldMode,
    /// A numerical decision fell between the versioned accept/reject bands.
    NumericalDecisionGrayZone,
    /// The checked peak-memory plan exceeded the supported capacity.
    CapacityExceeded,
    /// A backend candidate violated its backend-standard-form contract.
    BackendContractViolation,
    /// Backend-standard form passed, but physical recovery verification failed.
    RecoveryVerificationFailure,
    /// The convex feasible set is empty under an independently validated certificate.
    InfeasibleProblem,
    /// The convex objective is unbounded below along an independently validated recession ray.
    UnboundedProblem,
    /// Independently validated attempts reached mutually incompatible conclusions.
    NumericalConsistencyFailure,
    /// Numerical execution failed without proving a stronger diagnosis.
    NumericalFailure,
}

/// Stable preflight evidence for one unresolved Axial Normal.
#[derive(Debug, Clone, PartialEq)]
pub struct UnresolvedAxialNormalEvidence {
    source_id: SourceId,
    input_axis: Vector3,
    backend_invoked: bool,
}

impl UnresolvedAxialNormalEvidence {
    pub(crate) fn new(source_id: SourceId, input_axis: Vector3) -> Self {
        Self {
            source_id,
            input_axis,
            backend_invoked: false,
        }
    }

    /// Returns the unresolved caller-owned Axial Normal identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the normalized orientation retained from the original input axis.
    pub fn input_axis(&self) -> Vector3 {
        self.input_axis
    }

    /// Confirms that semantic preflight stopped before backend execution.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }
}

/// Stable proof that every supplied relation is invariant to a global constant shift.
#[derive(Debug, Clone, PartialEq)]
pub struct UnidentifiedAdditiveGaugeEvidence {
    source_ids: Box<[SourceId]>,
    group_ids: Box<[GroupId]>,
    backend_invoked: bool,
}

/// Stable canonical source association retained by certificate and recovery evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalEvidenceSource {
    source_id: SourceId,
    group_ids: Box<[GroupId]>,
    semantic_role: SemanticRolePath,
}

impl CanonicalEvidenceSource {
    pub(crate) fn new(
        source_id: SourceId,
        group_ids: Vec<GroupId>,
        semantic_role: SemanticRolePath,
    ) -> Self {
        Self {
            source_id,
            group_ids: group_ids.into(),
            semantic_role,
        }
    }

    /// Returns the caller-owned source identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns every stable semantic group associated with this evidence edge.
    pub fn group_ids(&self) -> &[GroupId] {
        &self.group_ids
    }

    /// Returns the canonical semantic role of this evidence edge.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }
}

/// One original canonical hard relation in a conflict witness.
#[derive(Debug, Clone, PartialEq)]
pub struct ConflictWitnessRelationEvidence {
    source: CanonicalEvidenceSource,
    multiplier: f64,
}

impl ConflictWitnessRelationEvidence {
    pub(crate) fn new(source: CanonicalEvidenceSource, multiplier: f64) -> Self {
        Self { source, multiplier }
    }

    /// Returns the original caller-owned relation provenance.
    pub fn source(&self) -> &CanonicalEvidenceSource {
        &self.source
    }

    /// Returns this relation's coefficient in the verified affine combination.
    pub fn multiplier(&self) -> f64 {
        self.multiplier
    }
}

/// A source-localized proof that a set of canonical hard relations is inconsistent.
#[derive(Debug, Clone, PartialEq)]
pub struct ConflictWitnessEvidence {
    relations: Box<[ConflictWitnessRelationEvidence]>,
    sources: Box<[CanonicalEvidenceSource]>,
    source_ids: Box<[SourceId]>,
    canonical_residual: f64,
    separation_margin: f64,
    residual_limit: f64,
    separation_limit: f64,
    provenance_verified: bool,
    backend_invoked: bool,
}

pub(crate) struct ConflictWitnessEvidenceParts {
    pub(crate) relations: Vec<ConflictWitnessRelationEvidence>,
    pub(crate) sources: Vec<CanonicalEvidenceSource>,
    pub(crate) canonical_residual: f64,
    pub(crate) separation_margin: f64,
    pub(crate) residual_limit: f64,
    pub(crate) separation_limit: f64,
    pub(crate) provenance_verified: bool,
    pub(crate) backend_invoked: bool,
}

impl ConflictWitnessEvidence {
    pub(crate) fn new(
        relations: Vec<ConflictWitnessRelationEvidence>,
        canonical_residual: f64,
        separation_margin: f64,
        residual_limit: f64,
        separation_limit: f64,
        provenance_verified: bool,
        backend_invoked: bool,
    ) -> Self {
        let sources = relations
            .iter()
            .map(|relation| relation.source.clone())
            .collect();
        Self::new_with_sources(ConflictWitnessEvidenceParts {
            relations,
            sources,
            canonical_residual,
            separation_margin,
            residual_limit,
            separation_limit,
            provenance_verified,
            backend_invoked,
        })
    }

    pub(crate) fn new_with_sources(parts: ConflictWitnessEvidenceParts) -> Self {
        let ConflictWitnessEvidenceParts {
            mut relations,
            mut sources,
            canonical_residual,
            separation_margin,
            residual_limit,
            separation_limit,
            provenance_verified,
            backend_invoked,
        } = parts;
        relations.sort_by(|left, right| {
            left.source
                .source_id
                .cmp(&right.source.source_id)
                .then_with(|| left.source.semantic_role.cmp(&right.source.semantic_role))
        });
        sources.sort_by(|left, right| {
            left.source_id
                .cmp(&right.source_id)
                .then_with(|| left.semantic_role.cmp(&right.semantic_role))
        });
        sources.dedup();
        let mut source_ids = sources
            .iter()
            .map(|source| source.source_id.clone())
            .collect::<Vec<_>>();
        source_ids.sort();
        source_ids.dedup();
        Self {
            relations: relations.into(),
            sources: sources.into(),
            source_ids: source_ids.into(),
            canonical_residual,
            separation_margin,
            residual_limit,
            separation_limit,
            provenance_verified,
            backend_invoked,
        }
    }

    /// Returns the original canonical relations and their proof multipliers.
    pub fn relations(&self) -> &[ConflictWitnessRelationEvidence] {
        &self.relations
    }

    /// Returns source provenance in stable SourceId order.
    pub fn sources(&self) -> &[CanonicalEvidenceSource] {
        &self.sources
    }

    /// Returns the distinct caller-owned SourceIds in stable order.
    pub fn source_ids(&self) -> &[SourceId] {
        &self.source_ids
    }

    /// Returns the maximum residual of the recomputed canonical combination.
    pub fn canonical_residual(&self) -> f64 {
        self.canonical_residual
    }

    /// Returns the strict separation proved by the recomputed targets.
    pub fn separation_margin(&self) -> f64 {
        self.separation_margin
    }

    /// Returns the canonical-combination residual acceptance limit.
    pub fn residual_limit(&self) -> f64 {
        self.residual_limit
    }

    /// Returns the minimum accepted strict-separation margin.
    pub fn separation_limit(&self) -> f64 {
        self.separation_limit
    }

    /// Reports whether every witness relation was recovered to original provenance.
    pub fn provenance_verified(&self) -> bool {
        self.provenance_verified
    }

    /// Reports whether a backend supplied the candidate proof coefficients.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }
}

/// Independently validated Farkas-ray evidence for convex infeasibility.
#[derive(Debug, Clone, PartialEq)]
pub struct InfeasibilityCertificateEvidence {
    finite: bool,
    normalized_ray_norm: f64,
    stationarity_residual: f64,
    dual_cone_violation: f64,
    separation_margin: f64,
    residual_limit: f64,
    separation_limit: f64,
    recovery_round_trip_error: f64,
    provenance_verified: bool,
    sources: Box<[CanonicalEvidenceSource]>,
    backend_invoked: bool,
}

impl InfeasibilityCertificateEvidence {
    pub(crate) fn new(parts: InfeasibilityCertificateEvidenceParts) -> Self {
        Self {
            finite: parts.finite,
            normalized_ray_norm: parts.normalized_ray_norm,
            stationarity_residual: parts.stationarity_residual,
            dual_cone_violation: parts.dual_cone_violation,
            separation_margin: parts.separation_margin,
            residual_limit: parts.residual_limit,
            separation_limit: parts.separation_limit,
            recovery_round_trip_error: parts.recovery_round_trip_error,
            provenance_verified: parts.provenance_verified,
            sources: parts.sources.into(),
            backend_invoked: parts.backend_invoked,
        }
    }

    /// Reports whether every retained certificate quantity is finite.
    pub fn finite(&self) -> bool {
        self.finite
    }
    /// Returns the infinity norm after deterministic ray normalization.
    pub fn normalized_ray_norm(&self) -> f64 {
        self.normalized_ray_norm
    }
    /// Returns the normalized `A^T z` residual.
    pub fn stationarity_residual(&self) -> f64 {
        self.stationarity_residual
    }
    /// Returns the largest violation of the dual cone.
    pub fn dual_cone_violation(&self) -> f64 {
        self.dual_cone_violation
    }
    /// Returns normalized strict separation `-b^T z`.
    pub fn separation_margin(&self) -> f64 {
        self.separation_margin
    }
    /// Returns the fixed residual and cone-violation limit.
    pub fn residual_limit(&self) -> f64 {
        self.residual_limit
    }
    /// Returns the fixed minimum strict-separation margin.
    pub fn separation_limit(&self) -> f64 {
        self.separation_limit
    }
    /// Returns the scaled-ray recovery round-trip error.
    pub fn recovery_round_trip_error(&self) -> f64 {
        self.recovery_round_trip_error
    }
    /// Reports whether the complete canonical/backend provenance map was verified.
    pub fn provenance_verified(&self) -> bool {
        self.provenance_verified
    }
    /// Returns stable source/group/role associations for the certificate proof.
    pub fn sources(&self) -> &[CanonicalEvidenceSource] {
        &self.sources
    }
    /// Reports that the backend supplied the candidate ray.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }
}

pub(crate) struct InfeasibilityCertificateEvidenceParts {
    pub(crate) finite: bool,
    pub(crate) normalized_ray_norm: f64,
    pub(crate) stationarity_residual: f64,
    pub(crate) dual_cone_violation: f64,
    pub(crate) separation_margin: f64,
    pub(crate) residual_limit: f64,
    pub(crate) separation_limit: f64,
    pub(crate) recovery_round_trip_error: f64,
    pub(crate) provenance_verified: bool,
    pub(crate) sources: Vec<CanonicalEvidenceSource>,
    pub(crate) backend_invoked: bool,
}

/// Independently recovered and validated recession-ray evidence for unboundedness.
#[derive(Debug, Clone, PartialEq)]
pub struct RecessionRayEvidence {
    finite: bool,
    normalized_ray_norm: f64,
    hessian_null_residual: f64,
    constraint_ray_violation: f64,
    descent_margin: f64,
    residual_limit: f64,
    separation_limit: f64,
    recovery_round_trip_error: f64,
    provenance_verified: bool,
    sources: Box<[CanonicalEvidenceSource]>,
    backend_invoked: bool,
}

impl RecessionRayEvidence {
    pub(crate) fn new(parts: RecessionRayEvidenceParts) -> Self {
        Self {
            finite: parts.finite,
            normalized_ray_norm: parts.normalized_ray_norm,
            hessian_null_residual: parts.hessian_null_residual,
            constraint_ray_violation: parts.constraint_ray_violation,
            descent_margin: parts.descent_margin,
            residual_limit: parts.residual_limit,
            separation_limit: parts.separation_limit,
            recovery_round_trip_error: parts.recovery_round_trip_error,
            provenance_verified: parts.provenance_verified,
            sources: parts.sources.into(),
            backend_invoked: parts.backend_invoked,
        }
    }

    /// Reports whether every retained recession quantity is finite.
    pub fn finite(&self) -> bool {
        self.finite
    }
    /// Returns the infinity norm after deterministic ray normalization.
    pub fn normalized_ray_norm(&self) -> f64 {
        self.normalized_ray_norm
    }
    /// Returns the normalized residual of the zero-curvature condition `P d = 0`.
    pub fn hessian_null_residual(&self) -> f64 {
        self.hessian_null_residual
    }
    /// Returns the largest normalized equality or recession-cone violation.
    pub fn constraint_ray_violation(&self) -> f64 {
        self.constraint_ray_violation
    }
    /// Returns normalized strict objective descent `-q^T d`.
    pub fn descent_margin(&self) -> f64 {
        self.descent_margin
    }
    /// Returns the fixed residual and ray-violation limit.
    pub fn residual_limit(&self) -> f64 {
        self.residual_limit
    }
    /// Returns the fixed minimum strict-descent margin.
    pub fn separation_limit(&self) -> f64 {
        self.separation_limit
    }
    /// Returns the scaled-ray recovery round-trip error.
    pub fn recovery_round_trip_error(&self) -> f64 {
        self.recovery_round_trip_error
    }
    /// Reports whether the complete canonical/backend provenance map was verified.
    pub fn provenance_verified(&self) -> bool {
        self.provenance_verified
    }
    /// Returns stable source/group/role associations for the recession proof.
    pub fn sources(&self) -> &[CanonicalEvidenceSource] {
        &self.sources
    }
    /// Reports that the backend supplied the candidate ray.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }
}

pub(crate) struct RecessionRayEvidenceParts {
    pub(crate) finite: bool,
    pub(crate) normalized_ray_norm: f64,
    pub(crate) hessian_null_residual: f64,
    pub(crate) constraint_ray_violation: f64,
    pub(crate) descent_margin: f64,
    pub(crate) residual_limit: f64,
    pub(crate) separation_limit: f64,
    pub(crate) recovery_round_trip_error: f64,
    pub(crate) provenance_verified: bool,
    pub(crate) sources: Vec<CanonicalEvidenceSource>,
    pub(crate) backend_invoked: bool,
}

impl UnidentifiedAdditiveGaugeEvidence {
    pub(crate) fn new(
        source_ids: Vec<SourceId>,
        group_ids: Vec<GroupId>,
        backend_invoked: bool,
    ) -> Self {
        Self {
            source_ids: source_ids.into(),
            group_ids: group_ids.into(),
            backend_invoked,
        }
    }

    /// Returns every shift-invariant caller source in stable order.
    pub fn source_ids(&self) -> &[SourceId] {
        &self.source_ids
    }

    /// Returns every affected semantic latent in stable group order.
    pub fn group_ids(&self) -> &[GroupId] {
        &self.group_ids
    }

    /// Reports whether a backend was called before this structural conclusion.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }
}

/// Stable proof that one shared-level latent has no informative relation.
#[derive(Debug, Clone, PartialEq)]
pub struct UninformativeSharedLevelSetEvidence {
    group_id: GroupId,
    member_source_id: SourceId,
    backend_invoked: bool,
}

impl UninformativeSharedLevelSetEvidence {
    pub(crate) fn new(
        group_id: GroupId,
        member_source_id: SourceId,
        backend_invoked: bool,
    ) -> Self {
        Self {
            group_id,
            member_source_id,
            backend_invoked,
        }
    }

    /// Returns the uninformative semantic latent's stable group identity.
    pub fn group_id(&self) -> &GroupId {
        &self.group_id
    }

    /// Returns the only member relation's caller-owned identity.
    pub fn member_source_id(&self) -> &SourceId {
        &self.member_source_id
    }

    /// Reports whether a backend was called before this structural conclusion.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }
}

/// Stable source and target evidence for a direct hard-input contradiction.
#[derive(Debug, Clone, PartialEq)]
pub struct DirectInputConflictEvidence {
    first_source: CanonicalEvidenceSource,
    second_source: CanonicalEvidenceSource,
    semantic_role: SemanticRolePath,
    first_target: f64,
    second_target: f64,
    canonical_residual: f64,
    separation_margin: f64,
    provenance_verified: bool,
}

impl DirectInputConflictEvidence {
    pub(crate) fn new_verified_same_lhs(
        first_source: CanonicalEvidenceSource,
        second_source: CanonicalEvidenceSource,
        first_target: f64,
        second_target: f64,
    ) -> Self {
        let semantic_role = second_source.semantic_role.clone();
        let raw_margin = (second_target - first_target).abs();
        let separation_margin = if raw_margin.is_finite() {
            raw_margin
        } else {
            let scale = first_target.abs().max(second_target.abs());
            (second_target / scale - first_target / scale).abs()
        };
        Self {
            first_source,
            second_source,
            semantic_role,
            first_target,
            second_target,
            canonical_residual: 0.0,
            separation_margin,
            provenance_verified: separation_margin > 0.0,
        }
    }

    /// Returns the first conflicting SourceId in stable order.
    pub fn first_source(&self) -> &SourceId {
        &self.first_source.source_id
    }

    /// Returns the second conflicting SourceId in stable order.
    pub fn second_source(&self) -> &SourceId {
        &self.second_source.source_id
    }

    /// Returns the conflicting scalar semantic component.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the first exact hard target.
    pub fn first_target(&self) -> f64 {
        self.first_target
    }

    /// Returns the incompatible exact hard target.
    pub fn second_target(&self) -> f64 {
        self.second_target
    }

    pub(crate) fn first_canonical_source(&self) -> &CanonicalEvidenceSource {
        &self.first_source
    }

    pub(crate) fn second_canonical_source(&self) -> &CanonicalEvidenceSource {
        &self.second_source
    }

    pub(crate) fn canonical_residual(&self) -> f64 {
        self.canonical_residual
    }

    pub(crate) fn separation_margin(&self) -> f64 {
        self.separation_margin
    }

    pub(crate) fn provenance_verified(&self) -> bool {
        self.provenance_verified
    }
}

/// Complete provenance for a contradiction proved by a hard-relation graph.
#[derive(Debug, Clone, PartialEq)]
pub struct RelationGraphConflictEvidence {
    proof_relations: Box<[ConflictWitnessRelationEvidence]>,
    source_ids: Box<[SourceId]>,
    group_ids: Box<[GroupId]>,
    semantic_role: SemanticRolePath,
    first_absolute_source: SourceId,
    first_absolute_target: f64,
    second_absolute_source: SourceId,
    second_absolute_target: f64,
    canonical_residual: f64,
    separation_margin: f64,
    provenance_verified: bool,
    backend_invoked: bool,
}

pub(crate) struct RelationGraphConflictEvidenceParts {
    pub(crate) proof_relations: Vec<ConflictWitnessRelationEvidence>,
    pub(crate) group_ids: Vec<GroupId>,
    pub(crate) semantic_role: SemanticRolePath,
    pub(crate) first_absolute_source: SourceId,
    pub(crate) first_absolute_target: f64,
    pub(crate) second_absolute_source: SourceId,
    pub(crate) second_absolute_target: f64,
    pub(crate) canonical_residual: f64,
    pub(crate) separation_margin: f64,
    pub(crate) provenance_verified: bool,
}

impl RelationGraphConflictEvidence {
    pub(crate) fn new(parts: RelationGraphConflictEvidenceParts) -> Self {
        let RelationGraphConflictEvidenceParts {
            proof_relations,
            group_ids,
            semantic_role,
            first_absolute_source,
            first_absolute_target,
            second_absolute_source,
            second_absolute_target,
            canonical_residual,
            separation_margin,
            provenance_verified,
        } = parts;
        let source_ids: Vec<_> = proof_relations
            .iter()
            .map(|relation| relation.source().source_id().clone())
            .collect();
        Self {
            proof_relations: proof_relations.into(),
            source_ids: source_ids.into(),
            group_ids: group_ids.into(),
            semantic_role,
            first_absolute_source,
            first_absolute_target,
            second_absolute_source,
            second_absolute_target,
            canonical_residual,
            separation_margin,
            provenance_verified,
            backend_invoked: false,
        }
    }

    pub(crate) fn proof_relations(&self) -> &[ConflictWitnessRelationEvidence] {
        &self.proof_relations
    }

    pub(crate) fn canonical_residual(&self) -> f64 {
        self.canonical_residual
    }

    pub(crate) fn separation_margin(&self) -> f64 {
        self.separation_margin
    }

    pub(crate) fn provenance_verified(&self) -> bool {
        self.provenance_verified
    }

    /// Returns every caller-owned source on the contradictory graph cycle.
    pub fn source_ids(&self) -> &[SourceId] {
        &self.source_ids
    }

    /// Returns every referenced semantic group on the contradictory cycle.
    pub fn group_ids(&self) -> &[GroupId] {
        &self.group_ids
    }

    /// Returns the semantic role of the relation that closed the cycle.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the source of the first incompatible absolute target.
    pub fn first_absolute_source(&self) -> &SourceId {
        &self.first_absolute_source
    }

    /// Returns the first incompatible absolute target without derived arithmetic.
    pub fn first_absolute_target(&self) -> f64 {
        self.first_absolute_target
    }

    /// Returns the source of the second incompatible absolute target.
    pub fn second_absolute_source(&self) -> &SourceId {
        &self.second_absolute_source
    }

    /// Returns the second incompatible absolute target without derived arithmetic.
    pub fn second_absolute_target(&self) -> f64 {
        self.second_absolute_target
    }

    /// Reports whether a backend was invoked before proving the contradiction.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }
}

/// One caller-owned source and its role in a shared-level-set conflict proof.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedLevelSetConflictSourceEvidence {
    source_id: SourceId,
    group_ids: Box<[GroupId]>,
    semantic_role: SemanticRolePath,
}

impl SharedLevelSetConflictSourceEvidence {
    pub(crate) fn new(
        source_id: SourceId,
        mut group_ids: Vec<GroupId>,
        semantic_role: SemanticRolePath,
    ) -> Self {
        group_ids.sort();
        group_ids.dedup();
        Self {
            source_id,
            group_ids: group_ids.into(),
            semantic_role,
        }
    }

    /// Returns the caller-owned source identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the shared level sets connected to this source in the proof.
    pub fn group_ids(&self) -> &[GroupId] {
        &self.group_ids
    }

    /// Returns this source's original semantic role.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }
}

/// Complete provenance for an impossible hard shared-level-set relation proof.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedLevelSetRelationConflictEvidence {
    source_provenance: Box<[SharedLevelSetConflictSourceEvidence]>,
    source_ids: Box<[SourceId]>,
    group_ids: Box<[GroupId]>,
    semantic_roles: Box<[SemanticRolePath]>,
    backend_invoked: bool,
    proof_multipliers: Option<Box<[f64]>>,
    canonical_residual: Option<f64>,
    separation_margin: Option<f64>,
    provenance_verified: bool,
}

impl SharedLevelSetRelationConflictEvidence {
    pub(crate) fn new(mut source_provenance: Vec<SharedLevelSetConflictSourceEvidence>) -> Self {
        source_provenance.sort_by(|left, right| left.source_id.cmp(&right.source_id));
        source_provenance.dedup_by(|left, right| left.source_id == right.source_id);
        let source_ids = source_provenance
            .iter()
            .map(|source| source.source_id.clone())
            .collect::<Vec<_>>();
        let semantic_roles = source_provenance
            .iter()
            .map(|source| source.semantic_role.clone())
            .collect::<Vec<_>>();
        let mut group_ids = source_provenance
            .iter()
            .flat_map(|source| source.group_ids.iter().cloned())
            .collect::<Vec<_>>();
        group_ids.sort();
        group_ids.dedup();
        Self {
            source_provenance: source_provenance.into(),
            source_ids: source_ids.into(),
            group_ids: group_ids.into(),
            semantic_roles: semantic_roles.into(),
            backend_invoked: false,
            proof_multipliers: None,
            canonical_residual: None,
            separation_margin: None,
            provenance_verified: false,
        }
    }

    pub(crate) fn new_with_canonical_witness(
        source_provenance: Vec<SharedLevelSetConflictSourceEvidence>,
        separation_margin: f64,
    ) -> Self {
        Self::new_with_verified_source_multipliers(
            source_provenance
                .into_iter()
                .map(|source| (source, 1.0))
                .collect(),
            separation_margin,
            0.0,
        )
    }

    pub(crate) fn new_with_verified_source_multipliers(
        mut source_provenance: Vec<(SharedLevelSetConflictSourceEvidence, f64)>,
        separation_margin: f64,
        canonical_residual: f64,
    ) -> Self {
        source_provenance.sort_by(|left, right| left.0.source_id.cmp(&right.0.source_id));
        let (sources, multipliers): (Vec<_>, Vec<_>) = source_provenance.into_iter().unzip();
        let mut evidence = Self::new(sources);
        evidence.provenance_verified = !multipliers.is_empty()
            && multipliers.iter().all(|multiplier| multiplier.is_finite())
            && canonical_residual == 0.0
            && separation_margin.is_finite()
            && separation_margin > 0.0;
        evidence.proof_multipliers = Some(multipliers.into());
        evidence.canonical_residual = Some(canonical_residual);
        evidence.separation_margin = Some(separation_margin);
        evidence
    }

    pub(crate) fn proof_multipliers(&self) -> Option<&[f64]> {
        self.proof_multipliers.as_deref()
    }

    pub(crate) fn canonical_residual(&self) -> Option<f64> {
        self.canonical_residual
    }

    pub(crate) fn separation_margin(&self) -> Option<f64> {
        self.separation_margin
    }

    pub(crate) fn provenance_verified(&self) -> bool {
        self.provenance_verified
    }

    /// Returns source/group/role associations for the complete proof.
    pub fn source_provenance(&self) -> &[SharedLevelSetConflictSourceEvidence] {
        &self.source_provenance
    }

    /// Returns every caller-owned relation source in stable order.
    pub fn source_ids(&self) -> &[SourceId] {
        &self.source_ids
    }

    /// Returns every shared level set participating in the proof in stable order.
    pub fn group_ids(&self) -> &[GroupId] {
        &self.group_ids
    }

    /// Returns the stable semantic role paired with each reported source.
    pub fn semantic_roles(&self) -> &[SemanticRolePath] {
        &self.semantic_roles
    }

    /// Reports whether a backend was invoked before proving the contradiction.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }
}

/// Physical dimension of one scalar hard-relation residual.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ResidualDimension {
    /// Scalar field-value units.
    FieldValue,
    /// Scalar field-value-per-length units.
    FieldValuePerLength,
}

/// Evidence for the implicit complete-Pi1 quotient construction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CubicQuotientConstructionEvidence {
    quotient_dimension: usize,
    householder_reflector_count: usize,
    congruence_pass_count: usize,
    householder_orthogonality_error: f64,
    canonical_response_round_trip_error: f64,
}

pub(crate) struct CubicQuotientConstructionEvidenceParts {
    pub(crate) quotient_dimension: usize,
    pub(crate) householder_reflector_count: usize,
    pub(crate) congruence_pass_count: usize,
    pub(crate) householder_orthogonality_error: f64,
    pub(crate) canonical_response_round_trip_error: f64,
}

impl CubicQuotientConstructionEvidence {
    pub(crate) fn new(parts: CubicQuotientConstructionEvidenceParts) -> Self {
        Self {
            quotient_dimension: parts.quotient_dimension,
            householder_reflector_count: parts.householder_reflector_count,
            congruence_pass_count: parts.congruence_pass_count,
            householder_orthogonality_error: parts.householder_orthogonality_error,
            canonical_response_round_trip_error: parts.canonical_response_round_trip_error,
        }
    }

    /// Returns the quotient dimension produced by the implicit construction.
    pub fn quotient_dimension(self) -> usize {
        self.quotient_dimension
    }

    /// Returns the number of Householder reflectors in the complete Pi1 QR.
    pub fn householder_reflector_count(self) -> usize {
        self.householder_reflector_count
    }

    /// Returns the number of full-matrix Householder applications used by the
    /// quotient congruence.
    pub fn congruence_pass_count(self) -> usize {
        self.congruence_pass_count
    }

    /// Returns the observed orthogonality defect of the complete Pi1
    /// Householder image.
    pub fn householder_orthogonality_error(self) -> f64 {
        self.householder_orthogonality_error
    }

    /// Returns the round-trip error for a canonical response through the
    /// implicit quotient coordinates.
    pub fn canonical_response_round_trip_error(self) -> f64 {
        self.canonical_response_round_trip_error
    }
}

/// Outward-rounded certificate interval for one unregularized quotient LLT pivot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CubicLltPivotInterval {
    lower: f64,
    upper: f64,
}

pub(crate) struct CubicLltPivotIntervalParts {
    pub(crate) lower: f64,
    pub(crate) upper: f64,
}

impl CubicLltPivotInterval {
    pub(crate) fn new(parts: CubicLltPivotIntervalParts) -> Self {
        Self {
            lower: parts.lower,
            upper: parts.upper,
        }
    }

    /// Returns the outward-rounded lower bound for this pivot.
    pub fn lower_bound(self) -> f64 {
        self.lower
    }

    /// Returns the outward-rounded upper bound for this pivot.
    pub fn upper_bound(self) -> f64 {
        self.upper
    }
}

/// Evidence for the verified, unregularized quotient LLT and reversible energy basis.
#[derive(Debug, Clone, PartialEq)]
pub struct CubicQuotientFactorizationEvidence {
    quotient_dimension: usize,
    retained_modes: usize,
    truncated_modes: usize,
    unregularized_llt_count: usize,
    full_spectrum_analysis_count: usize,
    normalized_backward_error: f64,
    pivot_intervals: Box<[CubicLltPivotInterval]>,
    field_energy_identity_error: f64,
    side_condition_error: f64,
    recovery_round_trip_error: f64,
    canonical_response_round_trip_error: f64,
    kernel_ridge_applied: bool,
    gram_jitter_applied: bool,
    mode_truncation_applied: bool,
    precision_rescue: Option<CubicPrecisionRescueEvidence>,
}

pub(crate) struct CubicQuotientFactorizationEvidenceParts {
    pub(crate) quotient_dimension: usize,
    pub(crate) retained_modes: usize,
    pub(crate) truncated_modes: usize,
    pub(crate) unregularized_llt_count: usize,
    pub(crate) full_spectrum_analysis_count: usize,
    pub(crate) normalized_backward_error: f64,
    pub(crate) pivot_intervals: Vec<CubicLltPivotInterval>,
    pub(crate) field_energy_identity_error: f64,
    pub(crate) side_condition_error: f64,
    pub(crate) recovery_round_trip_error: f64,
    pub(crate) canonical_response_round_trip_error: f64,
    pub(crate) kernel_ridge_applied: bool,
    pub(crate) gram_jitter_applied: bool,
    pub(crate) mode_truncation_applied: bool,
    pub(crate) precision_rescue: Option<CubicPrecisionRescueEvidence>,
}

impl CubicQuotientFactorizationEvidence {
    pub(crate) fn new(parts: CubicQuotientFactorizationEvidenceParts) -> Self {
        Self {
            quotient_dimension: parts.quotient_dimension,
            retained_modes: parts.retained_modes,
            truncated_modes: parts.truncated_modes,
            unregularized_llt_count: parts.unregularized_llt_count,
            full_spectrum_analysis_count: parts.full_spectrum_analysis_count,
            normalized_backward_error: parts.normalized_backward_error,
            pivot_intervals: parts.pivot_intervals.into(),
            field_energy_identity_error: parts.field_energy_identity_error,
            side_condition_error: parts.side_condition_error,
            recovery_round_trip_error: parts.recovery_round_trip_error,
            canonical_response_round_trip_error: parts.canonical_response_round_trip_error,
            kernel_ridge_applied: parts.kernel_ridge_applied,
            gram_jitter_applied: parts.gram_jitter_applied,
            mode_truncation_applied: parts.mode_truncation_applied,
            precision_rescue: parts.precision_rescue,
        }
    }

    /// Returns the full quotient dimension submitted to the LLT.
    pub fn quotient_dimension(&self) -> usize {
        self.quotient_dimension
    }

    /// Returns the number of positive effective field modes retained.
    pub fn retained_modes(&self) -> usize {
        self.retained_modes
    }

    /// Returns the number of quotient modes truncated from the representation.
    pub fn truncated_modes(&self) -> usize {
        self.truncated_modes
    }

    /// Returns the number of unregularized LLT calls on the complete quotient Gram matrix.
    pub fn unregularized_llt_count(&self) -> usize {
        self.unregularized_llt_count
    }

    /// Returns the number of full-spectrum analyses run after a successful LLT.
    pub fn full_spectrum_analysis_count(&self) -> usize {
        self.full_spectrum_analysis_count
    }

    /// Returns the scale-aware LLT backward certificate `eta_G`.
    pub fn normalized_backward_error(&self) -> f64 {
        self.normalized_backward_error
    }

    /// Returns every outward-rounded quotient pivot interval in factorization order.
    pub fn pivot_intervals(&self) -> &[CubicLltPivotInterval] {
        &self.pivot_intervals
    }

    /// Returns the recovered FieldEnergy-to-Euclidean identity defect.
    pub fn field_energy_identity_error(&self) -> f64 {
        self.field_energy_identity_error
    }

    /// Returns the complete-Pi1 side-condition defect after the energy change of basis.
    pub fn side_condition_error(&self) -> f64 {
        self.side_condition_error
    }

    /// Returns the solver-coordinate recovery round-trip defect.
    pub fn recovery_round_trip_error(&self) -> f64 {
        self.recovery_round_trip_error
    }

    /// Returns the canonical-response round-trip defect through energy coordinates.
    pub fn canonical_response_round_trip_error(&self) -> f64 {
        self.canonical_response_round_trip_error
    }

    /// Reports whether a kernel ridge altered the canonical pairing.
    pub fn kernel_ridge_applied(&self) -> bool {
        self.kernel_ridge_applied
    }

    /// Reports whether jitter altered the quotient Gram matrix.
    pub fn gram_jitter_applied(&self) -> bool {
        self.gram_jitter_applied
    }

    /// Reports whether any quotient mode was truncated.
    pub fn mode_truncation_applied(&self) -> bool {
        self.mode_truncation_applied
    }

    /// Returns the bounded upgrade performed for an ambiguous Schur block.
    pub fn precision_rescue(&self) -> Option<CubicPrecisionRescueEvidence> {
        self.precision_rescue
    }
}

/// Outcome of a bounded approximately 106-bit representation rescue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum CubicPrecisionRescueConclusion {
    /// Every upgraded mode was proved strictly positive.
    Positive,
    /// A reconstructible canonical algebraic zero was proved.
    AlgebraicZero,
    /// The upgraded Schur block contained negative curvature.
    NegativeCurvature,
    /// The upgraded interval still could not be separated from zero.
    NumericalDecisionGrayZone,
    /// More than the policy's 64-mode bound required upgrading.
    CapacityExceeded,
}

impl From<crate::precision_rescue::PrecisionRescueConclusion> for CubicPrecisionRescueConclusion {
    fn from(conclusion: crate::precision_rescue::PrecisionRescueConclusion) -> Self {
        match conclusion {
            crate::precision_rescue::PrecisionRescueConclusion::Positive => Self::Positive,
            crate::precision_rescue::PrecisionRescueConclusion::AlgebraicZero => {
                Self::AlgebraicZero
            }
            crate::precision_rescue::PrecisionRescueConclusion::NegativeCurvature => {
                Self::NegativeCurvature
            }
            crate::precision_rescue::PrecisionRescueConclusion::GrayZone => {
                Self::NumericalDecisionGrayZone
            }
            crate::precision_rescue::PrecisionRescueConclusion::CapacityExceeded => {
                Self::CapacityExceeded
            }
        }
    }
}

/// Auditable range, precision, and conclusion of one bounded rescue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CubicPrecisionRescueEvidence {
    first_mode: usize,
    mode_count: usize,
    precision_bits: u32,
    conclusion: CubicPrecisionRescueConclusion,
}

impl CubicPrecisionRescueEvidence {
    pub(crate) fn new(
        first_mode: usize,
        mode_count: usize,
        precision_bits: u32,
        conclusion: CubicPrecisionRescueConclusion,
    ) -> Self {
        Self {
            first_mode,
            mode_count,
            precision_bits,
            conclusion,
        }
    }

    /// Returns the first factorization-order mode in the upgraded block.
    pub fn first_mode(self) -> usize {
        self.first_mode
    }

    /// Returns the complete number of modes submitted to rescue.
    pub fn mode_count(self) -> usize {
        self.mode_count
    }

    /// Returns the arithmetic precision used by the bounded rescue.
    pub fn precision_bits(self) -> u32 {
        self.precision_bits
    }

    /// Returns the semantic conclusion reached after upgrading.
    pub fn conclusion(self) -> CubicPrecisionRescueConclusion {
        self.conclusion
    }
}

/// Complete Cubic representation analysis retained by a successful fit.
#[derive(Debug, Clone, PartialEq)]
pub struct CubicAnalysisEvidence {
    fitting_functional_count: usize,
    polynomial_dimension: usize,
    polynomial_rank: usize,
    quotient_construction: CubicQuotientConstructionEvidence,
    quotient_factorization: CubicQuotientFactorizationEvidence,
    polynomial_singular_values: Box<[f64]>,
    polynomial_rrqr_ratio: f64,
    polynomial_svd_ratio: f64,
    polynomial_rank_reject_ratio: f64,
    polynomial_rank_accept_ratio: f64,
    polynomial_precision_rescue: Option<CubicPrecisionRescueEvidence>,
    null_space_defect: f64,
    reduced_symmetry_defect: f64,
    reduced_symmetry_defect_limit: f64,
    reduced_largest_singular_value: f64,
    reduced_smallest_singular_value: f64,
    affine_reproduction_error: f64,
    solve_coordinate_length: f64,
    degenerate_extent: bool,
}

pub(crate) struct CubicAnalysisEvidenceParts {
    pub(crate) fitting_functional_count: usize,
    pub(crate) polynomial_dimension: usize,
    pub(crate) polynomial_rank: usize,
    pub(crate) quotient_construction: CubicQuotientConstructionEvidence,
    pub(crate) quotient_factorization: CubicQuotientFactorizationEvidence,
    pub(crate) polynomial_singular_values: Vec<f64>,
    pub(crate) polynomial_rrqr_ratio: f64,
    pub(crate) polynomial_svd_ratio: f64,
    pub(crate) polynomial_rank_reject_ratio: f64,
    pub(crate) polynomial_rank_accept_ratio: f64,
    pub(crate) polynomial_precision_rescue: Option<CubicPrecisionRescueEvidence>,
    pub(crate) null_space_defect: f64,
    pub(crate) reduced_symmetry_defect: f64,
    pub(crate) reduced_symmetry_defect_limit: f64,
    pub(crate) reduced_largest_singular_value: f64,
    pub(crate) reduced_smallest_singular_value: f64,
    pub(crate) affine_reproduction_error: f64,
    pub(crate) solve_coordinate_length: f64,
    pub(crate) degenerate_extent: bool,
}

impl CubicAnalysisEvidence {
    pub(crate) fn new(parts: CubicAnalysisEvidenceParts) -> Self {
        Self {
            fitting_functional_count: parts.fitting_functional_count,
            polynomial_dimension: parts.polynomial_dimension,
            polynomial_rank: parts.polynomial_rank,
            quotient_construction: parts.quotient_construction,
            quotient_factorization: parts.quotient_factorization,
            polynomial_singular_values: parts.polynomial_singular_values.into(),
            polynomial_rrqr_ratio: parts.polynomial_rrqr_ratio,
            polynomial_svd_ratio: parts.polynomial_svd_ratio,
            polynomial_rank_reject_ratio: parts.polynomial_rank_reject_ratio,
            polynomial_rank_accept_ratio: parts.polynomial_rank_accept_ratio,
            polynomial_precision_rescue: parts.polynomial_precision_rescue,
            null_space_defect: parts.null_space_defect,
            reduced_symmetry_defect: parts.reduced_symmetry_defect,
            reduced_symmetry_defect_limit: parts.reduced_symmetry_defect_limit,
            reduced_largest_singular_value: parts.reduced_largest_singular_value,
            reduced_smallest_singular_value: parts.reduced_smallest_singular_value,
            affine_reproduction_error: parts.affine_reproduction_error,
            solve_coordinate_length: parts.solve_coordinate_length,
            degenerate_extent: parts.degenerate_extent,
        }
    }

    /// Returns the number of unique fitting functionals in the representer span.
    pub fn fitting_functional_count(&self) -> usize {
        self.fitting_functional_count
    }

    /// Returns the complete Cubic polynomial-space dimension.
    pub fn polynomial_dimension(&self) -> usize {
        self.polynomial_dimension
    }

    /// Returns the accepted numerical rank of the polynomial pairing.
    pub fn polynomial_rank(&self) -> usize {
        self.polynomial_rank
    }

    /// Returns evidence for the implicit complete-Pi1 quotient construction.
    pub fn quotient_construction(&self) -> CubicQuotientConstructionEvidence {
        self.quotient_construction
    }

    /// Returns evidence for the verified energy-orthonormal quotient basis.
    pub fn quotient_factorization(&self) -> &CubicQuotientFactorizationEvidence {
        &self.quotient_factorization
    }

    /// Returns the polynomial-pairing singular values.
    pub fn polynomial_singular_values(&self) -> &[f64] {
        &self.polynomial_singular_values
    }

    /// Returns the polynomial RRQR rank ratio.
    pub fn polynomial_rrqr_ratio(&self) -> f64 {
        self.polynomial_rrqr_ratio
    }

    /// Returns the polynomial SVD rank ratio.
    pub fn polynomial_svd_ratio(&self) -> f64 {
        self.polynomial_svd_ratio
    }

    /// Returns the rank-rejection ratio boundary.
    pub fn polynomial_rank_reject_ratio(&self) -> f64 {
        self.polynomial_rank_reject_ratio
    }

    /// Returns the rank-acceptance ratio boundary.
    pub fn polynomial_rank_accept_ratio(&self) -> f64 {
        self.polynomial_rank_accept_ratio
    }

    /// Returns the bounded upgrade used to classify the complete Pi1 pairing.
    pub fn polynomial_precision_rescue(&self) -> Option<CubicPrecisionRescueEvidence> {
        self.polynomial_precision_rescue
    }

    /// Returns the null-space reconstruction defect.
    pub fn null_space_defect(&self) -> f64 {
        self.null_space_defect
    }

    /// Returns the reduced-pairing symmetry defect.
    pub fn reduced_symmetry_defect(&self) -> f64 {
        self.reduced_symmetry_defect
    }

    /// Returns the reduced-pairing symmetry-defect limit.
    pub fn reduced_symmetry_defect_limit(&self) -> f64 {
        self.reduced_symmetry_defect_limit
    }

    /// Returns a non-normative fixed-iteration estimate of the reduced
    /// pairing's largest singular value.
    pub fn reduced_largest_singular_value(&self) -> f64 {
        self.reduced_largest_singular_value
    }

    /// Returns a non-normative fixed-iteration estimate of the reduced
    /// pairing's smallest singular value.
    pub fn reduced_smallest_singular_value(&self) -> f64 {
        self.reduced_smallest_singular_value
    }

    /// Returns a non-normative largest-to-smallest risk estimate for the
    /// reduced pairing. This value never decides canonical rank.
    pub fn reduced_condition_estimate(&self) -> Option<f64> {
        let estimate = self.reduced_largest_singular_value / self.reduced_smallest_singular_value;
        estimate.is_finite().then_some(estimate)
    }

    /// Returns the complete-affine reproduction error.
    pub fn affine_reproduction_error(&self) -> f64 {
        self.affine_reproduction_error
    }

    /// Returns the characteristic solve-coordinate length.
    pub fn solve_coordinate_length(&self) -> f64 {
        self.solve_coordinate_length
    }

    /// Reports whether every support had zero geometric extent.
    pub fn degenerate_extent(&self) -> bool {
        self.degenerate_extent
    }
}

/// Algebraic object whose numerical rank was assessed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RankEvidenceDomain {
    /// Complete Cubic polynomial pairing.
    CubicPolynomialPairing,
    /// Cubic reduced kernel pairing.
    CubicReducedPairing,
    /// Symmetric augmented backend KKT matrix.
    BackendKkt,
}

/// Canonical field concept recovered from a numerical rank deficiency.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RankDeficiencyConcept {
    /// A mode in Cubic's complete affine polynomial space was not identified.
    CubicPi1FieldMode,
    /// A canonical Cubic quotient field mode had algebraically zero energy.
    CubicQuotientFieldMode,
}

/// Proof that a rank loss was interpreted in canonical field semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct InterpretableRankDeficiencyEvidence {
    concept: RankDeficiencyConcept,
    domain: RankEvidenceDomain,
    source_ids: Box<[SourceId]>,
    semantic_roles: Box<[SemanticRolePath]>,
    canonical_mode_residual: f64,
    canonical_mode_verified: bool,
    backend_invoked: bool,
    hidden_regularization_applied: bool,
}

pub(crate) struct InterpretableRankDeficiencyEvidenceParts {
    pub(crate) concept: RankDeficiencyConcept,
    pub(crate) domain: RankEvidenceDomain,
    pub(crate) source_ids: Vec<SourceId>,
    pub(crate) semantic_roles: Vec<SemanticRolePath>,
    pub(crate) canonical_mode_residual: f64,
    pub(crate) canonical_mode_verified: bool,
    pub(crate) backend_invoked: bool,
    pub(crate) hidden_regularization_applied: bool,
}

impl InterpretableRankDeficiencyEvidence {
    pub(crate) fn new(parts: InterpretableRankDeficiencyEvidenceParts) -> Self {
        Self {
            concept: parts.concept,
            domain: parts.domain,
            source_ids: parts.source_ids.into(),
            semantic_roles: parts.semantic_roles.into(),
            canonical_mode_residual: parts.canonical_mode_residual,
            canonical_mode_verified: parts.canonical_mode_verified,
            backend_invoked: parts.backend_invoked,
            hidden_regularization_applied: parts.hidden_regularization_applied,
        }
    }

    /// Returns the canonical field concept recovered from the rank evidence.
    pub fn concept(&self) -> RankDeficiencyConcept {
        self.concept
    }

    /// Returns the algebraic object whose null mode was interpreted.
    pub fn domain(&self) -> RankEvidenceDomain {
        self.domain
    }

    /// Returns every caller source participating in the interpreted mode.
    pub fn source_ids(&self) -> &[SourceId] {
        &self.source_ids
    }

    /// Returns the corresponding canonical semantic roles.
    pub fn semantic_roles(&self) -> &[SemanticRolePath] {
        &self.semantic_roles
    }

    /// Returns the independently recomputed canonical null-mode residual.
    pub fn canonical_mode_residual(&self) -> f64 {
        self.canonical_mode_residual
    }

    /// Reports whether the null mode was mapped and checked in canonical semantics.
    pub fn canonical_mode_verified(&self) -> bool {
        self.canonical_mode_verified
    }

    /// Reports whether a candidate-producing backend was invoked first.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }

    /// Reports whether any forbidden repair was used to alter the conclusion.
    pub fn hidden_regularization_applied(&self) -> bool {
        self.hidden_regularization_applied
    }
}

/// Result of a versioned rank decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RankDecision {
    /// The evidence lies above the acceptance boundary.
    FullRank,
    /// The evidence lies below the rejection boundary.
    RankDeficient,
    /// The evidence lies between the two decision boundaries.
    NumericalDecisionGrayZone,
}

/// Quantified rank evidence for a representation or backend matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct RankEvidence {
    domain: RankEvidenceDomain,
    dimension: usize,
    rank: Option<usize>,
    exact_zero_index: Option<usize>,
    rrqr_ratio: Option<f64>,
    singular_values: Box<[f64]>,
    svd_ratio: Option<f64>,
    reject_ratio: Option<f64>,
    accept_ratio: Option<f64>,
    decision: RankDecision,
    backend_invoked: bool,
}

pub(crate) struct RankEvidenceParts {
    pub(crate) domain: RankEvidenceDomain,
    pub(crate) dimension: usize,
    pub(crate) rank: Option<usize>,
    pub(crate) exact_zero_index: Option<usize>,
    pub(crate) rrqr_ratio: Option<f64>,
    pub(crate) singular_values: Vec<f64>,
    pub(crate) svd_ratio: Option<f64>,
    pub(crate) reject_ratio: Option<f64>,
    pub(crate) accept_ratio: Option<f64>,
    pub(crate) decision: RankDecision,
    pub(crate) backend_invoked: bool,
}

impl RankEvidence {
    pub(crate) fn new(parts: RankEvidenceParts) -> Self {
        Self {
            domain: parts.domain,
            dimension: parts.dimension,
            rank: parts.rank,
            exact_zero_index: parts.exact_zero_index,
            rrqr_ratio: parts.rrqr_ratio,
            singular_values: parts.singular_values.into(),
            svd_ratio: parts.svd_ratio,
            reject_ratio: parts.reject_ratio,
            accept_ratio: parts.accept_ratio,
            decision: parts.decision,
            backend_invoked: parts.backend_invoked,
        }
    }

    /// Returns the object whose rank was assessed.
    pub fn domain(&self) -> RankEvidenceDomain {
        self.domain
    }

    /// Returns the assessed square dimension.
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the estimated rank when the analysis produced one.
    pub fn rank(&self) -> Option<usize> {
        self.rank
    }

    /// Returns an exactly zero structural index, when found.
    pub fn exact_zero_index(&self) -> Option<usize> {
        self.exact_zero_index
    }

    /// Returns the RRQR decision ratio when available.
    pub fn rrqr_ratio(&self) -> Option<f64> {
        self.rrqr_ratio
    }

    /// Returns retained singular values in backend order.
    pub fn singular_values(&self) -> &[f64] {
        &self.singular_values
    }

    /// Returns the SVD decision ratio when available.
    pub fn svd_ratio(&self) -> Option<f64> {
        self.svd_ratio
    }

    /// Returns the rejection boundary when available.
    pub fn reject_ratio(&self) -> Option<f64> {
        self.reject_ratio
    }

    /// Returns the acceptance boundary when available.
    pub fn accept_ratio(&self) -> Option<f64> {
        self.accept_ratio
    }

    /// Returns the policy decision.
    pub fn decision(&self) -> RankDecision {
        self.decision
    }

    /// Reports whether the evidence proves full rank.
    pub fn is_full_rank(&self) -> bool {
        self.decision == RankDecision::FullRank
    }

    /// Reports whether a numerical backend had been invoked at this boundary.
    pub fn backend_invoked(&self) -> bool {
        self.backend_invoked
    }

    /// Returns a largest-to-smallest singular-value condition estimate.
    pub fn condition_estimate(&self) -> Option<f64> {
        let largest = self
            .singular_values
            .iter()
            .copied()
            .map(f64::abs)
            .reduce(f64::max)?;
        let smallest = self
            .singular_values
            .iter()
            .copied()
            .map(f64::abs)
            .reduce(f64::min)?;
        if smallest == 0.0 {
            return None;
        }
        let estimate = largest / smallest;
        estimate.is_finite().then_some(estimate)
    }
}

/// Positive, negative, and zero inertia counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InertiaCounts {
    positive: usize,
    negative: usize,
    zero: usize,
}

impl InertiaCounts {
    pub(crate) fn new(positive: usize, negative: usize, zero: usize) -> Self {
        Self {
            positive,
            negative,
            zero,
        }
    }

    /// Returns the positive eigenvalue count.
    pub fn positive(self) -> usize {
        self.positive
    }

    /// Returns the negative eigenvalue count.
    pub fn negative(self) -> usize {
        self.negative
    }

    /// Returns the zero eigenvalue count.
    pub fn zero(self) -> usize {
        self.zero
    }
}

/// Expected and observed inertia for the symmetric Equality KKT form.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InertiaEvidence {
    expected: InertiaCounts,
    observed: InertiaCounts,
    backend_invoked: bool,
}

impl InertiaEvidence {
    pub(crate) fn new(
        expected: InertiaCounts,
        observed: InertiaCounts,
        backend_invoked: bool,
    ) -> Self {
        Self {
            expected,
            observed,
            backend_invoked,
        }
    }

    /// Returns the convex Equality KKT inertia required by policy.
    pub fn expected(self) -> InertiaCounts {
        self.expected
    }

    /// Returns the independently observed inertia.
    pub fn observed(self) -> InertiaCounts {
        self.observed
    }

    /// Reports whether a candidate-producing backend was invoked.
    pub fn backend_invoked(self) -> bool {
        self.backend_invoked
    }
}

/// Physical Cubic side-condition assessment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SideConditionEvidence {
    components: [f64; 4],
    tolerances: [f64; 4],
    round_trip_error: f64,
}

impl SideConditionEvidence {
    pub(crate) fn new(components: [f64; 4], tolerances: [f64; 4], round_trip_error: f64) -> Self {
        Self {
            components,
            tolerances,
            round_trip_error,
        }
    }

    /// Returns constant/x/y/z side-condition components in physical coordinates.
    pub fn components(self) -> [f64; 4] {
        self.components
    }

    /// Returns the corresponding physical acceptance tolerances.
    pub fn tolerances(self) -> [f64; 4] {
        self.tolerances
    }

    /// Returns the forward/inverse side-condition round-trip error.
    pub fn round_trip_error(self) -> f64 {
        self.round_trip_error
    }
}

/// Independently counted participation and recovery evidence for every source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AllSourceRecoveryEvidence {
    canonical_hard_relation_count: usize,
    canonical_soft_relation_count: usize,
    participating_sources: Box<[SourceId]>,
    recovered_sources: Box<[SourceId]>,
    representer_count: usize,
    solver_relation_row_count: usize,
    recovery_edge_count: usize,
    verified: bool,
}

pub(crate) struct AllSourceRecoveryEvidenceParts {
    pub(crate) canonical_hard_relation_count: usize,
    pub(crate) canonical_soft_relation_count: usize,
    pub(crate) participating_sources: Vec<SourceId>,
    pub(crate) recovered_sources: Vec<SourceId>,
    pub(crate) representer_count: usize,
    pub(crate) solver_relation_row_count: usize,
    pub(crate) recovery_edge_count: usize,
    pub(crate) verified: bool,
}

impl AllSourceRecoveryEvidence {
    pub(crate) fn new(parts: AllSourceRecoveryEvidenceParts) -> Self {
        Self {
            canonical_hard_relation_count: parts.canonical_hard_relation_count,
            canonical_soft_relation_count: parts.canonical_soft_relation_count,
            participating_sources: parts.participating_sources.into(),
            recovered_sources: parts.recovered_sources.into(),
            representer_count: parts.representer_count,
            solver_relation_row_count: parts.solver_relation_row_count,
            recovery_edge_count: parts.recovery_edge_count,
            verified: parts.verified,
        }
    }

    /// Returns the independently counted canonical hard relations.
    pub fn canonical_hard_relation_count(&self) -> usize {
        self.canonical_hard_relation_count
    }

    /// Returns the independently counted canonical soft residual channels.
    pub fn canonical_soft_relation_count(&self) -> usize {
        self.canonical_soft_relation_count
    }

    /// Returns every SourceId whose canonical relation participates in the problem.
    pub fn participating_sources(&self) -> &[SourceId] {
        &self.participating_sources
    }

    /// Returns every SourceId reached through a verified recovery edge.
    pub fn recovered_sources(&self) -> &[SourceId] {
        &self.recovered_sources
    }

    /// Returns the independently constructed physical representer count.
    pub fn representer_count(&self) -> usize {
        self.representer_count
    }

    /// Returns relation rows retained by the solver-independent form.
    pub fn solver_relation_row_count(&self) -> usize {
        self.solver_relation_row_count
    }

    /// Returns source-bearing canonical-to-recovery edges.
    pub fn recovery_edge_count(&self) -> usize {
        self.recovery_edge_count
    }

    /// Reports complete source coverage and objective association verification.
    pub fn verified(&self) -> bool {
        self.verified
    }
}

/// Complete physical Recover-and-Verify acceptance evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct CanonicalAcceptanceEvidence {
    accepted: bool,
    backend_standard_form_verified: bool,
    recovery_finite: bool,
    provenance_verified: bool,
    side_condition: Option<SideConditionEvidence>,
    field_value_hard_residual_max: Option<f64>,
    field_value_per_length_hard_residual_max: Option<f64>,
    polynomial_round_trip_error: Option<f64>,
    field_coefficient_round_trip_error: Option<f64>,
    field_energy_round_trip_error: Option<f64>,
    whitening_round_trip_error: Option<f64>,
    objective_round_trip_error: Option<f64>,
    objective_verified: bool,
    tolerance_round_trip_error: Option<f64>,
    hard_affine_inequality_violation_max: Option<f64>,
    backend_standard_form_residual: Option<f64>,
    physical_convex_residual: Option<ConvexResidualEvidence>,
    scaling_round_trip_error: Option<f64>,
    reduction_round_trip_error: Option<f64>,
    backend_internal_scaling_round_trip_error: Option<f64>,
}

pub(crate) struct CanonicalAcceptanceEvidenceParts {
    pub(crate) accepted: bool,
    pub(crate) backend_standard_form_verified: bool,
    pub(crate) recovery_finite: bool,
    pub(crate) provenance_verified: bool,
    pub(crate) side_condition: Option<SideConditionEvidence>,
    pub(crate) hard_residual_maxima: Option<(f64, f64)>,
    pub(crate) polynomial_round_trip_error: Option<f64>,
    pub(crate) field_coefficient_round_trip_error: Option<f64>,
    pub(crate) field_energy_round_trip_error: Option<f64>,
    pub(crate) whitening_round_trip_error: Option<f64>,
    pub(crate) objective_round_trip_error: Option<f64>,
    pub(crate) objective_verified: bool,
    pub(crate) tolerance_round_trip_error: Option<f64>,
    pub(crate) hard_affine_inequality_violation_max: Option<f64>,
    pub(crate) backend_standard_form_residual: Option<f64>,
    pub(crate) physical_convex_residual: Option<ConvexResidualEvidence>,
    pub(crate) scaling_round_trip_error: Option<f64>,
    pub(crate) reduction_round_trip_error: Option<f64>,
    pub(crate) backend_internal_scaling_round_trip_error: Option<f64>,
}

impl CanonicalAcceptanceEvidence {
    pub(crate) fn new(parts: CanonicalAcceptanceEvidenceParts) -> Self {
        Self {
            accepted: parts.accepted,
            backend_standard_form_verified: parts.backend_standard_form_verified,
            recovery_finite: parts.recovery_finite,
            provenance_verified: parts.provenance_verified,
            side_condition: parts.side_condition,
            field_value_hard_residual_max: parts.hard_residual_maxima.map(|maxima| maxima.0),
            field_value_per_length_hard_residual_max: parts
                .hard_residual_maxima
                .map(|maxima| maxima.1),
            polynomial_round_trip_error: parts.polynomial_round_trip_error,
            field_coefficient_round_trip_error: parts.field_coefficient_round_trip_error,
            field_energy_round_trip_error: parts.field_energy_round_trip_error,
            whitening_round_trip_error: parts.whitening_round_trip_error,
            objective_round_trip_error: parts.objective_round_trip_error,
            objective_verified: parts.objective_verified,
            tolerance_round_trip_error: parts.tolerance_round_trip_error,
            hard_affine_inequality_violation_max: parts.hard_affine_inequality_violation_max,
            backend_standard_form_residual: parts.backend_standard_form_residual,
            physical_convex_residual: parts.physical_convex_residual,
            scaling_round_trip_error: parts.scaling_round_trip_error,
            reduction_round_trip_error: parts.reduction_round_trip_error,
            backend_internal_scaling_round_trip_error: parts
                .backend_internal_scaling_round_trip_error,
        }
    }

    /// Reports whether every canonical acceptance check passed.
    pub fn accepted(&self) -> bool {
        self.accepted
    }

    /// Reports whether the recovered backend-standard-form contract passed.
    pub fn backend_standard_form_verified(&self) -> bool {
        self.backend_standard_form_verified
    }

    /// Reports whether every recovered physical quantity was finite.
    pub fn recovery_finite(&self) -> bool {
        self.recovery_finite
    }

    /// Reports whether canonical provenance survived assembly and recovery.
    pub fn provenance_verified(&self) -> bool {
        self.provenance_verified
    }

    /// Returns the complete Cubic side-condition assessment when reached.
    pub fn side_condition(&self) -> Option<SideConditionEvidence> {
        self.side_condition
    }

    /// Returns the maximum hard field-value residual when reached.
    pub fn field_value_hard_residual_max(&self) -> Option<f64> {
        self.field_value_hard_residual_max
    }

    /// Returns the maximum hard derivative residual when reached.
    pub fn field_value_per_length_hard_residual_max(&self) -> Option<f64> {
        self.field_value_per_length_hard_residual_max
    }

    /// Returns polynomial recovery round-trip error when reached.
    pub fn polynomial_round_trip_error(&self) -> Option<f64> {
        self.polynomial_round_trip_error
    }

    /// Returns field-coefficient recovery round-trip error when reached.
    pub fn field_coefficient_round_trip_error(&self) -> Option<f64> {
        self.field_coefficient_round_trip_error
    }

    /// Returns FieldEnergy recovery round-trip error when reached.
    pub fn field_energy_round_trip_error(&self) -> Option<f64> {
        self.field_energy_round_trip_error
    }

    /// Returns the maximum whitening recovery round-trip error when reached.
    pub fn whitening_round_trip_error(&self) -> Option<f64> {
        self.whitening_round_trip_error
    }

    /// Returns physical/standard objective recovery round-trip error when
    /// reached.
    pub fn objective_round_trip_error(&self) -> Option<f64> {
        self.objective_round_trip_error
    }

    /// Reports whether the independently recomputed physical objective passed
    /// the recovery contract.
    pub fn objective_verified(&self) -> bool {
        self.objective_verified
    }

    /// Returns relation-tolerance recovery round-trip error when reached.
    pub fn tolerance_round_trip_error(&self) -> Option<f64> {
        self.tolerance_round_trip_error
    }

    /// Returns the largest recovered hard affine-inequality violation.
    pub fn hard_affine_inequality_violation_max(&self) -> Option<f64> {
        self.hard_affine_inequality_violation_max
    }

    /// Returns the recovered backend-standard-form residual for a QP candidate.
    pub fn backend_standard_form_residual(&self) -> Option<f64> {
        self.backend_standard_form_residual
    }

    /// Returns the five-part convex residual envelope recomputed after
    /// recovering the candidate into the physical canonical QP coordinates.
    pub fn physical_convex_residual(&self) -> Option<ConvexResidualEvidence> {
        self.physical_convex_residual
    }

    /// Returns the QP candidate scaling recovery round-trip error.
    pub fn scaling_round_trip_error(&self) -> Option<f64> {
        self.scaling_round_trip_error
    }

    /// Returns the QP null-space reduction round-trip error.
    pub fn reduction_round_trip_error(&self) -> Option<f64> {
        self.reduction_round_trip_error
    }

    /// Returns Clarabel's independently checked internal-scaling round-trip error.
    pub fn backend_internal_scaling_round_trip_error(&self) -> Option<f64> {
        self.backend_internal_scaling_round_trip_error
    }
}

/// Kind of checked capacity failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum CapacityFailureKind {
    /// Checked size arithmetic overflowed.
    ArithmeticOverflow,
    /// The planned peak exceeded the fixed capacity limit.
    LimitExceeded,
}

/// Capacity evidence retained when allocation planning rejects a fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CapacityEvidence {
    kind: CapacityFailureKind,
    limit_bytes: u64,
    planned_peak_bytes: Option<u64>,
    large_allocation_attempted: bool,
    backend_invocation_attempted: bool,
}

impl CapacityEvidence {
    pub(crate) fn new(
        kind: CapacityFailureKind,
        limit_bytes: u64,
        planned_peak_bytes: Option<u64>,
        large_allocation_attempted: bool,
        backend_invocation_attempted: bool,
    ) -> Self {
        Self {
            kind,
            limit_bytes,
            planned_peak_bytes,
            large_allocation_attempted,
            backend_invocation_attempted,
        }
    }

    /// Returns why the checked capacity plan failed.
    pub fn kind(self) -> CapacityFailureKind {
        self.kind
    }

    /// Returns the fixed peak-memory limit.
    pub fn limit_bytes(self) -> u64 {
        self.limit_bytes
    }

    /// Returns the planned peak when it was representable.
    pub fn planned_peak_bytes(self) -> Option<u64> {
        self.planned_peak_bytes
    }

    /// Reports whether a large allocation had been attempted.
    pub fn large_allocation_attempted(self) -> bool {
        self.large_allocation_attempted
    }

    /// Reports whether a numerical backend had been invoked.
    pub fn backend_invocation_attempted(self) -> bool {
        self.backend_invocation_attempted
    }
}

/// Failure point within Cubic representation or Equality KKT analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AnalysisFailureStage {
    /// Complete Cubic polynomial-rank analysis.
    CubicPolynomialRank,
    /// Cholesky positivity check of the reduced Cubic pairing.
    CubicReducedCholesky,
    /// Inertia check of the reduced Cubic pairing.
    CubicReducedInertia,
    /// Spectral analysis of the reduced Cubic pairing.
    CubicReducedSpectrum,
    /// Backend KKT rank confirmation.
    BackendRankConfirmation,
    /// Backend KKT inertia analysis.
    BackendInertia,
    /// Backend rank-analysis workspace.
    BackendRankWorkspace,
    /// Backend inertia-analysis workspace.
    BackendInertiaWorkspace,
    /// Primary symmetric-indefinite factorization workspace.
    BackendFactorWorkspace,
    /// Primary linear-solve workspace.
    BackendSolveWorkspace,
    /// Full-SVD rescue workspace.
    BackendSvdRescueWorkspace,
}

/// Reason a finite input could not define an invertible solve-coordinate map.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SolveCoordinateFailureReason {
    /// The bounding-box center was not finite.
    BoundingBoxCenterNotFinite,
    /// The characteristic support length was not finite.
    CharacteristicLengthNotFinite,
    /// Cubing the characteristic length did not produce an invertible scale.
    FieldRecoveryScaleNotInvertible,
    /// Transforming a canonical functional produced a non-finite value.
    StandardFunctionalNotFinite,
}

/// Physical or representation contract quantity that exceeded its limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AnalysisContractQuantity {
    /// Null-space reconstruction defect.
    NullSpaceDefect,
    /// Householder orthogonality defect.
    HouseholderOrthogonalityError,
    /// Canonical response quotient-coordinate round-trip error.
    CanonicalResponseRoundTripError,
    /// Unregularized quotient LLT backward certificate.
    QuotientLltBackwardError,
    /// FieldEnergy identity defect in solver-facing quotient coordinates.
    QuotientFieldEnergyIdentityError,
    /// Complete-Pi1 side-condition defect after energy orthonormalization.
    QuotientSideConditionError,
    /// Solver-to-canonical basis recovery round-trip error.
    QuotientRecoveryRoundTripError,
    /// Canonical response round-trip error through the energy basis.
    QuotientBasisResponseRoundTripError,
    /// Reduced-pairing symmetry defect.
    ReducedSymmetryDefect,
    /// Complete-affine reproduction error.
    AffineReproductionError,
}

/// Backend-standard-form input array.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum BackendInputField {
    /// Primal Hessian.
    Hessian,
    /// Equality Jacobian.
    EqualityJacobian,
    /// Stationarity right-hand side.
    StationarityRightHandSide,
    /// Equality right-hand side.
    EqualityRightHandSide,
}

/// Deterministic KKT scaling failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ScalingFailureReason {
    /// Matrix, right-hand-side, and dimension shapes disagreed.
    InvalidShape,
    /// A row had exactly zero norm.
    ZeroNorm { index: usize },
    /// A row norm was non-finite.
    NonFiniteNorm { index: usize },
}

/// Structured evidence for a numerical-analysis failure without a candidate.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum AnalysisFailureEvidence {
    /// The representer span was empty after canonical construction.
    EmptyRepresenterSpan,
    /// A checked physical-to-solve coordinate transform was not invertible.
    InvalidSolveCoordinateTransform {
        reason: SolveCoordinateFailureReason,
        backend_invoked: bool,
    },
    /// A numerical analysis failed without a stronger rank or inertia result.
    NumericalAnalysis {
        stage: AnalysisFailureStage,
        backend_invoked: bool,
    },
    /// A checked workspace allocation failed.
    WorkspaceAllocation {
        stage: AnalysisFailureStage,
        bytes: u64,
        alignment: usize,
        backend_invoked: bool,
    },
    /// Null-space construction workspace could not be allocated.
    NullSpaceWorkspaceAllocation,
    /// Complete-Pi1 classification remained unresolved after bounded rescue.
    PolynomialPrecisionRescue {
        rescue: CubicPrecisionRescueEvidence,
        backend_invoked: bool,
    },
    /// An f64 quotient LLT pivot requires the bounded precision-rescue stage
    /// before it can be classified.
    QuotientPivotRequiresPrecisionRescue {
        quotient_dimension: usize,
        pivot_index: usize,
        interval: Option<CubicLltPivotInterval>,
        backend_invoked: bool,
    },
    /// Bounded double-double rescue completed without a positive or algebraic-zero proof.
    QuotientPrecisionRescue {
        quotient_dimension: usize,
        rescue: CubicPrecisionRescueEvidence,
        backend_invoked: bool,
    },
    /// A quotient LLT pivot was reliably non-positive in f64. This is a
    /// numerical representation failure, not proof of canonical rank loss.
    QuotientFactorizationNotPositive {
        quotient_dimension: usize,
        pivot_index: usize,
        interval: CubicLltPivotInterval,
        backend_invoked: bool,
    },
    /// A quantified representation contract exceeded its policy limit.
    ContractThresholdExceeded {
        quantity: AnalysisContractQuantity,
        observed: f64,
        limit: f64,
    },
    /// A backend-standard-form array had the wrong checked length.
    InvalidBackendInputLength {
        field: BackendInputField,
        expected: usize,
        actual: usize,
    },
    /// A backend-standard-form input contained a non-finite value.
    NonFiniteBackendInput {
        field: BackendInputField,
        index: usize,
    },
    /// Deterministic diagonal KKT scaling could not be formed.
    ScalingFailure { reason: ScalingFailureReason },
}

/// Backend termination evidence, distinct from [`ProblemDiagnosis`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SolveAttemptTermination {
    /// The backend produced a nominal-accuracy candidate.
    CandidateProduced,
    /// The backend produced a reduced-accuracy candidate.
    ReducedAccuracyCandidateProduced,
    /// The backend produced evidence suggesting primal infeasibility.
    PrimalInfeasibilityCandidate,
    /// The backend produced evidence suggesting dual infeasibility.
    DualInfeasibilityCandidate,
    /// The configured iteration or resource limit was reached.
    LimitReached,
    /// The backend stopped after insufficient progress.
    InsufficientProgress,
    /// A configured callback requested termination.
    CallbackTermination,
    /// The backend attempt stopped on a numerical error.
    NumericalError,
}

/// Auditable identity of the backend used by one fit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BackendFingerprint {
    schema_version: u32,
    crate_name: Box<str>,
    crate_version: Box<str>,
    legacy_feature_slots: [Box<str>; 2],
    enabled_features: Box<[Box<str>]>,
    algorithm: Box<str>,
    target_arch: Box<str>,
    target_os: Box<str>,
    requested_threads: usize,
    actual_threads: usize,
}

pub(crate) struct BackendFingerprintParts {
    pub(crate) schema_version: u32,
    pub(crate) crate_name: &'static str,
    pub(crate) crate_version: &'static str,
    pub(crate) features: Vec<&'static str>,
    pub(crate) algorithm: &'static str,
    pub(crate) target_arch: &'static str,
    pub(crate) target_os: &'static str,
    pub(crate) requested_threads: usize,
    pub(crate) actual_threads: usize,
}

impl BackendFingerprint {
    pub(crate) fn new(parts: BackendFingerprintParts) -> Self {
        let legacy_feature_slots = [
            parts.features.first().copied().unwrap_or_default().into(),
            parts.features.get(1).copied().unwrap_or_default().into(),
        ];
        Self {
            schema_version: parts.schema_version,
            crate_name: parts.crate_name.into(),
            crate_version: parts.crate_version.into(),
            legacy_feature_slots,
            enabled_features: parts.features.into_iter().map(Into::into).collect(),
            algorithm: parts.algorithm.into(),
            target_arch: parts.target_arch.into(),
            target_os: parts.target_os.into(),
            requested_threads: parts.requested_threads,
            actual_threads: parts.actual_threads,
        }
    }

    /// Returns the fingerprint schema version.
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Returns the backend crate name.
    pub fn crate_name(&self) -> &str {
        &self.crate_name
    }

    /// Returns the backend crate version.
    pub fn crate_version(&self) -> &str {
        &self.crate_version
    }

    /// Returns the two legacy backend feature slots.
    ///
    /// This fixed-width view is retained for v0.1 source compatibility. A
    /// backend with fewer than two enabled features leaves trailing slots
    /// empty; use [`Self::enabled_features`] for the exact feature set.
    pub fn features(&self) -> [&str; 2] {
        self.legacy_feature_slots
            .each_ref()
            .map(|feature| feature.as_ref())
    }

    /// Iterates over the exact enabled backend features.
    pub fn enabled_features(&self) -> impl ExactSizeIterator<Item = &str> + '_ {
        self.enabled_features.iter().map(|feature| feature.as_ref())
    }

    /// Returns the resolved backend algorithm.
    pub fn algorithm(&self) -> &str {
        &self.algorithm
    }

    /// Returns the target architecture recorded by the adapter.
    pub fn target_arch(&self) -> &str {
        &self.target_arch
    }

    /// Returns the target operating system recorded by the adapter.
    pub fn target_os(&self) -> &str {
        &self.target_os
    }

    /// Returns the requested thread count recorded by the backend adapter.
    pub fn requested_threads(&self) -> usize {
        self.requested_threads
    }

    /// Returns the actual thread count recorded by the backend adapter.
    pub fn actual_threads(&self) -> usize {
        self.actual_threads
    }
}

/// Algorithm selected for one bounded backend attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SolveAttemptKind {
    /// Symmetric-indefinite Bunch-Kaufman factorization with refinement.
    BunchKaufmanRefinement,
    /// Full-SVD rescue after the primary candidate was rejected.
    FullSvdRescue,
    /// Clarabel's standard deterministic QP profile.
    ClarabelStandard,
    /// Clarabel's robust deterministic QP retry profile.
    ClarabelRobust,
}

/// Resolved settings for one backend attempt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BackendAttemptSettings {
    kind: SolveAttemptKind,
    pivoting: Option<Box<str>>,
    block_size: Option<usize>,
    parallelism_threshold: Option<usize>,
    factor_workspace_source: Option<Box<str>>,
    maximum_refinement_steps: Option<usize>,
    settings_id: Option<Box<str>>,
    left_vectors: Option<Box<str>>,
    right_vectors: Option<Box<str>>,
}

impl BackendAttemptSettings {
    pub(crate) fn lblt(
        pivoting: impl Into<Box<str>>,
        block_size: usize,
        parallelism_threshold: usize,
        factor_workspace_source: impl Into<Box<str>>,
        maximum_refinement_steps: usize,
    ) -> Self {
        Self {
            kind: SolveAttemptKind::BunchKaufmanRefinement,
            pivoting: Some(pivoting.into()),
            block_size: Some(block_size),
            parallelism_threshold: Some(parallelism_threshold),
            factor_workspace_source: Some(factor_workspace_source.into()),
            maximum_refinement_steps: Some(maximum_refinement_steps),
            settings_id: None,
            left_vectors: None,
            right_vectors: None,
        }
    }

    pub(crate) fn full_svd(
        settings_id: impl Into<Box<str>>,
        left_vectors: impl Into<Box<str>>,
        right_vectors: impl Into<Box<str>>,
    ) -> Self {
        Self {
            kind: SolveAttemptKind::FullSvdRescue,
            pivoting: None,
            block_size: None,
            parallelism_threshold: None,
            factor_workspace_source: None,
            maximum_refinement_steps: None,
            settings_id: Some(settings_id.into()),
            left_vectors: Some(left_vectors.into()),
            right_vectors: Some(right_vectors.into()),
        }
    }

    pub(crate) fn clarabel(kind: SolveAttemptKind, settings_id: impl Into<Box<str>>) -> Self {
        debug_assert!(matches!(
            kind,
            SolveAttemptKind::ClarabelStandard | SolveAttemptKind::ClarabelRobust
        ));
        Self {
            kind,
            pivoting: None,
            block_size: None,
            parallelism_threshold: None,
            factor_workspace_source: None,
            maximum_refinement_steps: None,
            settings_id: Some(settings_id.into()),
            left_vectors: None,
            right_vectors: None,
        }
    }

    /// Returns the algorithm family these settings configure.
    pub fn kind(&self) -> SolveAttemptKind {
        self.kind
    }

    /// Returns the pivoting strategy for a factorization attempt.
    pub fn pivoting(&self) -> Option<&str> {
        self.pivoting.as_deref()
    }

    /// Returns the factorization block size when applicable.
    pub fn block_size(&self) -> Option<usize> {
        self.block_size
    }

    /// Returns the backend parallelism threshold when applicable.
    pub fn parallelism_threshold(&self) -> Option<usize> {
        self.parallelism_threshold
    }

    /// Returns how factorization workspace was provisioned.
    pub fn factor_workspace_source(&self) -> Option<&str> {
        self.factor_workspace_source.as_deref()
    }

    /// Returns the maximum refinement-step budget when applicable.
    pub fn maximum_refinement_steps(&self) -> Option<usize> {
        self.maximum_refinement_steps
    }

    /// Returns the backend settings identity for an SVD rescue.
    pub fn settings_id(&self) -> Option<&str> {
        self.settings_id.as_deref()
    }

    /// Returns the requested left-singular-vector mode when applicable.
    pub fn left_vectors(&self) -> Option<&str> {
        self.left_vectors.as_deref()
    }

    /// Returns the requested right-singular-vector mode when applicable.
    pub fn right_vectors(&self) -> Option<&str> {
        self.right_vectors.as_deref()
    }
}

/// Ruiz scaling summary applied to one backend attempt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScalingSummary {
    method: Box<str>,
    rounds: usize,
    saturated_outside_target: usize,
}

impl ScalingSummary {
    pub(crate) fn new(
        method: impl Into<Box<str>>,
        rounds: usize,
        saturated_outside_target: usize,
    ) -> Self {
        Self {
            method: method.into(),
            rounds,
            saturated_outside_target,
        }
    }

    /// Returns the stable scaling-method identity.
    pub fn method(&self) -> &str {
        &self.method
    }

    /// Returns the number of completed fixed scaling rounds.
    pub fn rounds(&self) -> usize {
        self.rounds
    }

    /// Returns how many rows saturated outside the target norm band.
    pub fn saturated_outside_target(&self) -> usize {
        self.saturated_outside_target
    }
}

/// Complete dimensionless residual evidence for a backend candidate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearResidualEvidence {
    infinity_norm: f64,
    matrix_infinity_norm: f64,
    solution_infinity_norm: f64,
    rhs_infinity_norm: f64,
    normalized_backward_error: f64,
}

impl LinearResidualEvidence {
    pub(crate) fn new(values: [f64; 5]) -> Self {
        Self {
            infinity_norm: values[0],
            matrix_infinity_norm: values[1],
            solution_infinity_norm: values[2],
            rhs_infinity_norm: values[3],
            normalized_backward_error: values[4],
        }
    }

    /// Returns the candidate residual infinity norm.
    pub fn infinity_norm(self) -> f64 {
        self.infinity_norm
    }

    /// Returns the backend-standard-form matrix infinity norm.
    pub fn matrix_infinity_norm(self) -> f64 {
        self.matrix_infinity_norm
    }

    /// Returns the backend candidate infinity norm.
    pub fn solution_infinity_norm(self) -> f64 {
        self.solution_infinity_norm
    }

    /// Returns the right-hand-side infinity norm.
    pub fn rhs_infinity_norm(self) -> f64 {
        self.rhs_infinity_norm
    }

    /// Returns the dimensionless normalized backward error.
    pub fn normalized_backward_error(self) -> f64 {
        self.normalized_backward_error
    }
}

/// Category of a rejected backend attempt or execution path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AttemptFailureCategory {
    /// The backend returned a non-finite candidate.
    NonFiniteCandidate,
    /// The candidate exceeded the normalized backward-error limit.
    BackwardErrorExceeded,
    /// Scaling recovery exceeded its round-trip limit.
    ScalingRoundTripExceeded,
    /// The backend decomposition itself failed numerically.
    BackendDecompositionFailure,
    /// A QP candidate exceeded a convex residual acceptance limit.
    ConvexResidualExceeded,
    /// The backend did not honor the one-thread contract.
    ThreadContractViolation,
    /// The backend identity or complete settings fingerprint changed.
    BackendFingerprintMismatch,
    /// The termination did not carry a verified candidate or certificate.
    UnverifiedTermination,
    /// A claimed primal-infeasibility ray failed independent validation.
    InvalidInfeasibilityCertificate,
    /// A claimed dual-infeasibility ray failed independent validation.
    InvalidRecessionCertificate,
    /// Independently validated attempts reached incompatible conclusions.
    InconsistentValidatedConclusions,
}

/// Independently recomputed dimensionless residuals for a convex candidate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConvexResidualEvidence {
    primal: f64,
    dual: f64,
    stationarity: f64,
    complementarity: f64,
    relative_gap: f64,
}

impl ConvexResidualEvidence {
    pub(crate) fn new(parts: ConvexResidualEvidenceParts) -> Self {
        Self {
            primal: parts.primal,
            dual: parts.dual,
            stationarity: parts.stationarity,
            complementarity: parts.complementarity,
            relative_gap: parts.relative_gap,
        }
    }

    /// Returns the scaled primal residual.
    pub fn primal(self) -> f64 {
        self.primal
    }
    /// Returns the scaled dual-cone residual.
    pub fn dual(self) -> f64 {
        self.dual
    }
    /// Returns the scaled stationarity residual.
    pub fn stationarity(self) -> f64 {
        self.stationarity
    }
    /// Returns the scaled complementarity residual.
    pub fn complementarity(self) -> f64 {
        self.complementarity
    }
    /// Returns the relative primal-dual gap.
    pub fn relative_gap(self) -> f64 {
        self.relative_gap
    }
}

pub(crate) struct ConvexResidualEvidenceParts {
    pub(crate) primal: f64,
    pub(crate) dual: f64,
    pub(crate) stationarity: f64,
    pub(crate) complementarity: f64,
    pub(crate) relative_gap: f64,
}

/// Structured reason why an attempt or execution path was rejected.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AttemptFailureEvidence {
    category: AttemptFailureCategory,
    observed: Option<f64>,
    limit: Option<f64>,
}

impl AttemptFailureEvidence {
    pub(crate) fn new(
        category: AttemptFailureCategory,
        observed: Option<f64>,
        limit: Option<f64>,
    ) -> Self {
        Self {
            category,
            observed,
            limit,
        }
    }

    /// Returns the stable failure category.
    pub fn category(self) -> AttemptFailureCategory {
        self.category
    }

    /// Returns the observed value for thresholded failures.
    pub fn observed(self) -> Option<f64> {
        self.observed
    }

    /// Returns the applied acceptance limit for thresholded failures.
    pub fn limit(self) -> Option<f64> {
        self.limit
    }
}

/// Evidence retained for one bounded backend attempt.
#[derive(Debug, Clone, PartialEq)]
pub struct SolveAttemptRecord {
    sequence: usize,
    kind: SolveAttemptKind,
    termination: SolveAttemptTermination,
    settings: BackendAttemptSettings,
    scaling: ScalingSummary,
    scaling_round_trip_error: Option<f64>,
    refinement_steps: usize,
    residual: Option<LinearResidualEvidence>,
    convex_residual: Option<ConvexResidualEvidence>,
    certificate_present: bool,
    failure_reason: Option<AttemptFailureEvidence>,
    backend_fingerprint: BackendFingerprint,
}

pub(crate) struct SolveAttemptRecordParts {
    pub(crate) sequence: usize,
    pub(crate) kind: SolveAttemptKind,
    pub(crate) termination: SolveAttemptTermination,
    pub(crate) settings: BackendAttemptSettings,
    pub(crate) scaling: ScalingSummary,
    pub(crate) scaling_round_trip_error: Option<f64>,
    pub(crate) refinement_steps: usize,
    pub(crate) residual: Option<LinearResidualEvidence>,
    pub(crate) convex_residual: Option<ConvexResidualEvidence>,
    pub(crate) certificate_present: bool,
    pub(crate) failure_reason: Option<AttemptFailureEvidence>,
    pub(crate) backend_fingerprint: BackendFingerprint,
}

impl SolveAttemptRecord {
    pub(crate) fn new(parts: SolveAttemptRecordParts) -> Self {
        Self {
            sequence: parts.sequence,
            kind: parts.kind,
            termination: parts.termination,
            settings: parts.settings,
            scaling: parts.scaling,
            scaling_round_trip_error: parts.scaling_round_trip_error,
            refinement_steps: parts.refinement_steps,
            residual: parts.residual,
            convex_residual: parts.convex_residual,
            certificate_present: parts.certificate_present,
            failure_reason: parts.failure_reason,
            backend_fingerprint: parts.backend_fingerprint,
        }
    }

    /// Returns the deterministic attempt sequence number.
    pub fn sequence(&self) -> usize {
        self.sequence
    }

    /// Returns the algorithm family used for this attempt.
    pub fn kind(&self) -> SolveAttemptKind {
        self.kind
    }

    /// Returns backend termination evidence for this attempt.
    pub fn termination(&self) -> SolveAttemptTermination {
        self.termination
    }

    /// Returns the complete resolved attempt settings.
    pub fn settings(&self) -> &BackendAttemptSettings {
        &self.settings
    }

    /// Returns the scaling summary applied to this attempt.
    pub fn scaling(&self) -> &ScalingSummary {
        &self.scaling
    }

    /// Returns the scaling forward/inverse round-trip error when retained.
    pub fn scaling_round_trip_error(&self) -> Option<f64> {
        self.scaling_round_trip_error
    }

    /// Returns the number of completed iterative-refinement steps.
    pub fn refinement_steps(&self) -> usize {
        self.refinement_steps
    }

    /// Returns complete residual evidence when a candidate was produced.
    pub fn residual(&self) -> Option<LinearResidualEvidence> {
        self.residual
    }

    /// Returns independently recomputed convex residuals for a QP attempt.
    pub fn convex_residual(&self) -> Option<ConvexResidualEvidence> {
        self.convex_residual
    }

    /// Returns normalized backward error when the attempt produced a candidate.
    pub fn normalized_backward_error(&self) -> Option<f64> {
        self.residual
            .map(|evidence| evidence.normalized_backward_error())
    }

    /// Reports whether this attempt retained a backend certificate.
    pub fn certificate_present(&self) -> bool {
        self.certificate_present
    }

    /// Returns structured rejection evidence for a failed attempt.
    pub fn failure_reason(&self) -> Option<AttemptFailureEvidence> {
        self.failure_reason
    }

    /// Returns the adapter-recorded backend identity and settings.
    pub fn backend_fingerprint(&self) -> &BackendFingerprint {
        &self.backend_fingerprint
    }
}

/// Structured physical evidence retained when recovery rejects a candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct RecoveryVerificationEvidence {
    reasons: Vec<RecoveryVerificationReason>,
    side_condition: Option<SideConditionEvidence>,
    field_value_hard_residual_max: Option<f64>,
    field_value_per_length_hard_residual_max: Option<f64>,
    polynomial_round_trip_error: Option<f64>,
    field_coefficient_round_trip_error: Option<f64>,
    field_energy_round_trip_error: Option<f64>,
    whitening_round_trip_error: Option<f64>,
    objective_round_trip_error: Option<f64>,
    tolerance_round_trip_error: Option<f64>,
    backend_standard_form_residual: Option<f64>,
    reduction_round_trip_error: Option<f64>,
    scaling_round_trip_error: Option<f64>,
    sources: Box<[CanonicalEvidenceSource]>,
    no_model_produced: bool,
}

pub(crate) struct RecoveryVerificationEvidenceParts {
    pub(crate) reasons: Vec<RecoveryVerificationReason>,
    pub(crate) side_condition: Option<SideConditionEvidence>,
    pub(crate) hard_residual_maxima: Option<(f64, f64)>,
    pub(crate) polynomial_round_trip_error: Option<f64>,
    pub(crate) field_coefficient_round_trip_error: Option<f64>,
    pub(crate) field_energy_round_trip_error: Option<f64>,
    pub(crate) whitening_round_trip_error: Option<f64>,
    pub(crate) objective_round_trip_error: Option<f64>,
    pub(crate) tolerance_round_trip_error: Option<f64>,
    pub(crate) backend_standard_form_residual: Option<f64>,
    pub(crate) reduction_round_trip_error: Option<f64>,
    pub(crate) scaling_round_trip_error: Option<f64>,
    pub(crate) sources: Vec<CanonicalEvidenceSource>,
    pub(crate) no_model_produced: bool,
}

impl RecoveryVerificationEvidence {
    pub(crate) fn new(parts: RecoveryVerificationEvidenceParts) -> Self {
        Self {
            reasons: parts.reasons,
            side_condition: parts.side_condition,
            field_value_hard_residual_max: parts.hard_residual_maxima.map(|maxima| maxima.0),
            field_value_per_length_hard_residual_max: parts
                .hard_residual_maxima
                .map(|maxima| maxima.1),
            polynomial_round_trip_error: parts.polynomial_round_trip_error,
            field_coefficient_round_trip_error: parts.field_coefficient_round_trip_error,
            field_energy_round_trip_error: parts.field_energy_round_trip_error,
            whitening_round_trip_error: parts.whitening_round_trip_error,
            objective_round_trip_error: parts.objective_round_trip_error,
            tolerance_round_trip_error: parts.tolerance_round_trip_error,
            backend_standard_form_residual: parts.backend_standard_form_residual,
            reduction_round_trip_error: parts.reduction_round_trip_error,
            scaling_round_trip_error: parts.scaling_round_trip_error,
            sources: parts.sources.into(),
            no_model_produced: parts.no_model_produced,
        }
    }

    /// Returns every recovery rejection reason in deterministic order.
    pub fn reasons(&self) -> &[RecoveryVerificationReason] {
        &self.reasons
    }

    /// Returns the complete physical side-condition assessment when reached.
    pub fn side_condition(&self) -> Option<SideConditionEvidence> {
        self.side_condition
    }

    /// Returns the maximum hard field-value residual when recovery reached it.
    pub fn field_value_hard_residual_max(&self) -> Option<f64> {
        self.field_value_hard_residual_max
    }

    /// Returns the maximum hard derivative residual when recovery reached it.
    pub fn field_value_per_length_hard_residual_max(&self) -> Option<f64> {
        self.field_value_per_length_hard_residual_max
    }

    /// Returns polynomial round-trip error when available.
    pub fn polynomial_round_trip_error(&self) -> Option<f64> {
        self.polynomial_round_trip_error
    }

    /// Returns field-coefficient round-trip error when available.
    pub fn field_coefficient_round_trip_error(&self) -> Option<f64> {
        self.field_coefficient_round_trip_error
    }

    /// Returns FieldEnergy round-trip error when available.
    pub fn field_energy_round_trip_error(&self) -> Option<f64> {
        self.field_energy_round_trip_error
    }

    /// Returns the maximum whitening recovery round-trip error when available.
    pub fn whitening_round_trip_error(&self) -> Option<f64> {
        self.whitening_round_trip_error
    }

    /// Returns physical/standard objective round-trip error when available.
    pub fn objective_round_trip_error(&self) -> Option<f64> {
        self.objective_round_trip_error
    }

    /// Returns relation-tolerance round-trip error when available.
    pub fn tolerance_round_trip_error(&self) -> Option<f64> {
        self.tolerance_round_trip_error
    }

    /// Returns the recovered backend-standard-form residual for a QP candidate.
    pub fn backend_standard_form_residual(&self) -> Option<f64> {
        self.backend_standard_form_residual
    }

    /// Returns the QP null-space reduction recovery error when available.
    pub fn reduction_round_trip_error(&self) -> Option<f64> {
        self.reduction_round_trip_error
    }

    /// Returns the GeoRBF QP scaling recovery error when available.
    pub fn scaling_round_trip_error(&self) -> Option<f64> {
        self.scaling_round_trip_error
    }

    /// Returns stable source/group/role associations reached during recovery.
    pub fn sources(&self) -> &[CanonicalEvidenceSource] {
        &self.sources
    }

    /// Confirms that rejection produced no public model.
    pub fn no_model_produced(&self) -> bool {
        self.no_model_produced
    }
}
