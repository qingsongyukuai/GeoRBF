//! Typed public fit diagnoses and backend-attempt evidence.

pub use crate::cubic_equality::RecoveryVerificationFailureReason as RecoveryVerificationReason;
pub use crate::functional::SemanticRolePath;
use crate::functional::{GroupId, SourceId};
pub use crate::numerical::NumericalPolicyId;

/// GeoRBF's semantic conclusion for a failed fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProblemDiagnosis {
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
    /// Numerical execution failed without proving a stronger diagnosis.
    NumericalFailure,
}

/// Stable proof that every supplied relation is invariant to a global constant shift.
#[derive(Debug, Clone, PartialEq)]
pub struct UnidentifiedAdditiveGaugeEvidence {
    source_ids: Box<[SourceId]>,
    group_ids: Box<[GroupId]>,
    backend_invoked: bool,
}

/// Independently validated Farkas-ray evidence for convex infeasibility.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InfeasibilityCertificateEvidence {
    finite: bool,
    normalized_ray_norm: f64,
    stationarity_residual: f64,
    dual_cone_violation: f64,
    separation_margin: f64,
    residual_limit: f64,
    separation_limit: f64,
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
            backend_invoked: parts.backend_invoked,
        }
    }

    /// Reports whether every retained certificate quantity is finite.
    pub fn finite(self) -> bool {
        self.finite
    }
    /// Returns the infinity norm after deterministic ray normalization.
    pub fn normalized_ray_norm(self) -> f64 {
        self.normalized_ray_norm
    }
    /// Returns the normalized `A^T z` residual.
    pub fn stationarity_residual(self) -> f64 {
        self.stationarity_residual
    }
    /// Returns the largest violation of the dual cone.
    pub fn dual_cone_violation(self) -> f64 {
        self.dual_cone_violation
    }
    /// Returns normalized strict separation `-b^T z`.
    pub fn separation_margin(self) -> f64 {
        self.separation_margin
    }
    /// Returns the fixed residual and cone-violation limit.
    pub fn residual_limit(self) -> f64 {
        self.residual_limit
    }
    /// Returns the fixed minimum strict-separation margin.
    pub fn separation_limit(self) -> f64 {
        self.separation_limit
    }
    /// Reports that the backend supplied the candidate ray.
    pub fn backend_invoked(self) -> bool {
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
    first_source: SourceId,
    second_source: SourceId,
    semantic_role: SemanticRolePath,
    first_target: f64,
    second_target: f64,
}

impl DirectInputConflictEvidence {
    pub(crate) fn new(
        first_source: SourceId,
        second_source: SourceId,
        semantic_role: SemanticRolePath,
        first_target: f64,
        second_target: f64,
    ) -> Self {
        Self {
            first_source,
            second_source,
            semantic_role,
            first_target,
            second_target,
        }
    }

    /// Returns the first conflicting SourceId in stable order.
    pub fn first_source(&self) -> &SourceId {
        &self.first_source
    }

    /// Returns the second conflicting SourceId in stable order.
    pub fn second_source(&self) -> &SourceId {
        &self.second_source
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
}

/// Complete provenance for a contradiction proved by a hard-relation graph.
#[derive(Debug, Clone, PartialEq)]
pub struct RelationGraphConflictEvidence {
    source_ids: Box<[SourceId]>,
    group_ids: Box<[GroupId]>,
    semantic_role: SemanticRolePath,
    first_absolute_source: SourceId,
    first_absolute_target: f64,
    second_absolute_source: SourceId,
    second_absolute_target: f64,
    backend_invoked: bool,
}

impl RelationGraphConflictEvidence {
    pub(crate) fn new(
        source_ids: Vec<SourceId>,
        group_ids: Vec<GroupId>,
        semantic_role: SemanticRolePath,
        first_absolute_source: SourceId,
        first_absolute_target: f64,
        second_absolute_source: SourceId,
        second_absolute_target: f64,
    ) -> Self {
        Self {
            source_ids: source_ids.into(),
            group_ids: group_ids.into(),
            semantic_role,
            first_absolute_source,
            first_absolute_target,
            second_absolute_source,
            second_absolute_target,
            backend_invoked: false,
        }
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

/// Physical dimension of one scalar hard-relation residual.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ResidualDimension {
    /// Scalar field-value units.
    FieldValue,
    /// Scalar field-value-per-length units.
    FieldValuePerLength,
}

/// Complete Cubic representation analysis retained by a successful fit.
#[derive(Debug, Clone, PartialEq)]
pub struct CubicAnalysisEvidence {
    fitting_functional_count: usize,
    polynomial_dimension: usize,
    polynomial_rank: usize,
    polynomial_singular_values: Box<[f64]>,
    polynomial_rrqr_ratio: f64,
    polynomial_svd_ratio: f64,
    polynomial_rank_reject_ratio: f64,
    polynomial_rank_accept_ratio: f64,
    null_space_defect: f64,
    reduced_symmetry_defect: f64,
    reduced_symmetry_defect_limit: f64,
    reduced_smallest_singular_value: f64,
    affine_reproduction_error: f64,
    solve_coordinate_length: f64,
    degenerate_extent: bool,
}

pub(crate) struct CubicAnalysisEvidenceParts {
    pub(crate) fitting_functional_count: usize,
    pub(crate) polynomial_dimension: usize,
    pub(crate) polynomial_rank: usize,
    pub(crate) polynomial_singular_values: Vec<f64>,
    pub(crate) polynomial_rrqr_ratio: f64,
    pub(crate) polynomial_svd_ratio: f64,
    pub(crate) polynomial_rank_reject_ratio: f64,
    pub(crate) polynomial_rank_accept_ratio: f64,
    pub(crate) null_space_defect: f64,
    pub(crate) reduced_symmetry_defect: f64,
    pub(crate) reduced_symmetry_defect_limit: f64,
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
            polynomial_singular_values: parts.polynomial_singular_values.into(),
            polynomial_rrqr_ratio: parts.polynomial_rrqr_ratio,
            polynomial_svd_ratio: parts.polynomial_svd_ratio,
            polynomial_rank_reject_ratio: parts.polynomial_rank_reject_ratio,
            polynomial_rank_accept_ratio: parts.polynomial_rank_accept_ratio,
            null_space_defect: parts.null_space_defect,
            reduced_symmetry_defect: parts.reduced_symmetry_defect,
            reduced_symmetry_defect_limit: parts.reduced_symmetry_defect_limit,
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

    /// Returns the reduced pairing's smallest singular value.
    pub fn reduced_smallest_singular_value(&self) -> f64 {
        self.reduced_smallest_singular_value
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

/// Complete physical Recover-and-Verify acceptance evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct CanonicalAcceptanceEvidence {
    accepted: bool,
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
}

pub(crate) struct CanonicalAcceptanceEvidenceParts {
    pub(crate) accepted: bool,
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
}

impl CanonicalAcceptanceEvidence {
    pub(crate) fn new(parts: CanonicalAcceptanceEvidenceParts) -> Self {
        Self {
            accepted: parts.accepted,
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
        }
    }

    /// Reports whether every canonical acceptance check passed.
    pub fn accepted(&self) -> bool {
        self.accepted
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
    features: Box<[Box<str>]>,
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
        Self {
            schema_version: parts.schema_version,
            crate_name: parts.crate_name.into(),
            crate_version: parts.crate_version.into(),
            features: parts.features.into_iter().map(Into::into).collect(),
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

    /// Returns the exact enabled backend features.
    pub fn features(&self) -> &[Box<str>] {
        &self.features
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

    /// Confirms that rejection produced no public model.
    pub fn no_model_produced(&self) -> bool {
        self.no_model_produced
    }
}
