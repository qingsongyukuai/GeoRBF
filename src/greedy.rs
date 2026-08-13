//! Frozen Greedy call-chain audit and deterministic public-fit evidence.
//!
//! Sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - `surfe_lib/surfe_api.cpp` (`SetGreedyAlgorithm`, `ComputeInterpolant`);
//! - `surfe_lib/modeling_methods.{h,cpp}` (`get_method`, `run_greedy_algorithm`);
//! - `surfe_lib/modelling_input.{h,cpp}` (four residual selectors);
//! - the five model `get_minimial_and_excluded_input`, `measure_residuals`, and
//!   `append_greedy_input` overrides.
//!
//! Frozen `SetGreedyAlgorithm` stores the request, but `ComputeInterpolant`
//! never reads `use_greedy` and never calls `run_greedy_algorithm`. GeoRBF
//! therefore records a zero-round trace for the public fit instead of exposing
//! the source-only loop, its TODO hooks, or its undefined tangent-selector
//! branch as supported behavior.

use crate::{ModelType, Parameters};

/// Source-level state of one frozen model hook.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum GreedyHookBody {
    /// A substantive source body exists, but has no `Surfe_API` call edge.
    Implemented,
    /// The override is an inline `return true` marked TODO (or equivalent).
    TodoStub,
}

/// Auditable frozen Greedy hook classification for one model kind.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct GreedyModelAudit {
    pub model: ModelType,
    /// Model constructed by `GRBF_Modelling_Methods::get_method` if the
    /// source-only loop is invoked directly.
    pub greedy_factory_model: ModelType,
    /// Always false: `Surfe_API::ComputeInterpolant` has no call edge to the
    /// Greedy loop for any model kind.
    pub reachable_from_surfe_api: bool,
    pub minimal: GreedyHookBody,
    pub residual: GreedyHookBody,
    pub append: GreedyHookBody,
}

/// Complete five-model source audit in [`ModelType::ALL`] order.
pub const GREEDY_MODEL_AUDIT: [GreedyModelAudit; 5] = [
    GreedyModelAudit {
        model: ModelType::SingleSurface,
        greedy_factory_model: ModelType::SingleSurface,
        reachable_from_surfe_api: false,
        minimal: GreedyHookBody::Implemented,
        residual: GreedyHookBody::Implemented,
        append: GreedyHookBody::Implemented,
    },
    GreedyModelAudit {
        model: ModelType::LajaunieApproach,
        greedy_factory_model: ModelType::LajaunieApproach,
        reachable_from_surfe_api: false,
        minimal: GreedyHookBody::Implemented,
        residual: GreedyHookBody::Implemented,
        append: GreedyHookBody::Implemented,
    },
    GreedyModelAudit {
        model: ModelType::StratigraphicHorizons,
        greedy_factory_model: ModelType::StratigraphicHorizons,
        reachable_from_surfe_api: false,
        minimal: GreedyHookBody::TodoStub,
        residual: GreedyHookBody::TodoStub,
        append: GreedyHookBody::TodoStub,
    },
    GreedyModelAudit {
        model: ModelType::ContinuousProperty,
        greedy_factory_model: ModelType::ContinuousProperty,
        reachable_from_surfe_api: false,
        minimal: GreedyHookBody::TodoStub,
        residual: GreedyHookBody::Implemented,
        append: GreedyHookBody::Implemented,
    },
    GreedyModelAudit {
        model: ModelType::VectorField,
        greedy_factory_model: ModelType::ContinuousProperty,
        reachable_from_surfe_api: false,
        minimal: GreedyHookBody::TodoStub,
        residual: GreedyHookBody::TodoStub,
        append: GreedyHookBody::TodoStub,
    },
];

/// Evidence schema for a Greedy iteration.
///
/// Public frozen fits never produce an instance. The fields make the absence
/// of selection and residual work explicit rather than representing it as an
/// undocumented empty optimization result.
#[derive(Clone, Debug, PartialEq)]
pub struct GreedyRoundEvidence {
    pub iteration: usize,
    pub selected_inequality_indices: Vec<usize>,
    pub selected_interface_indices: Vec<usize>,
    pub selected_planar_indices: Vec<usize>,
    pub selected_tangent_indices: Vec<usize>,
    pub inequality_residuals: Vec<bool>,
    pub interface_residuals: Vec<f64>,
    pub planar_residuals: Vec<f64>,
    pub tangent_residuals: Vec<f64>,
}

/// Why a public fitted-model Greedy trace stopped.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum GreedyStopReason {
    /// `Surfe_API::ComputeInterpolant` bypasses `run_greedy_algorithm`.
    NotCalledBySurfeApi,
}

/// Deterministic evidence for the Greedy request attached to a public fit.
#[derive(Clone, Debug, PartialEq)]
pub struct GreedyTrace {
    /// Stored `Parameters::use_greedy`. The frozen setter writes `true` even
    /// when its boolean argument is false.
    pub stored_use_greedy: bool,
    pub interface_uncertainty: f64,
    pub angular_uncertainty: f64,
    /// Always empty for the frozen public call chain.
    pub rounds: Vec<GreedyRoundEvidence>,
    pub stop_reason: GreedyStopReason,
}

impl GreedyTrace {
    pub(crate) fn public_fit(parameters: &Parameters) -> Self {
        Self {
            stored_use_greedy: parameters.use_greedy,
            interface_uncertainty: parameters.interface_uncertainty,
            angular_uncertainty: parameters.angular_uncertainty,
            rounds: Vec::new(),
            stop_reason: GreedyStopReason::NotCalledBySurfeApi,
        }
    }
}
