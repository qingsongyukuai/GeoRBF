use std::error::Error;

use georbf_cubic_cpd_recovery_spike::{
    CLARABEL_VERSION, FAER_VERSION, run_counterexamples, run_manufactured_experiment,
};

fn main() -> Result<(), Box<dyn Error>> {
    let evidence = run_manufactured_experiment()?;
    let counterexamples = run_counterexamples()?;

    println!("backend.faer={FAER_VERSION}");
    println!("backend.clarabel={CLARABEL_VERSION}");
    println!("cpd.functionals={}", evidence.cpd.functional_count);
    println!(
        "cpd.polynomial.dimension={}",
        evidence.cpd.polynomial_dimension
    );
    println!("cpd.polynomial.rank={}", evidence.cpd.polynomial_rank);
    println!(
        "cpd.polynomial.singular_values={:?}",
        evidence.cpd.singular_values
    );
    println!(
        "cpd.null_space.defect={:.17e}",
        evidence.cpd.null_space_defect
    );
    println!(
        "cpd.reduced.symmetry_defect={:.17e}",
        evidence.cpd.reduced_symmetry_defect
    );
    println!(
        "cpd.reduced.symmetry_limit={:.17e}",
        evidence.cpd.symmetry_defect_limit
    );
    println!(
        "cpd.reduced.smallest_eigenvalue={:.17e}",
        evidence.cpd.reduced_smallest_eigenvalue
    );
    println!(
        "cpd.affine_reproduction_error={:.17e}",
        evidence.cpd.affine_reproduction_error
    );
    println!(
        "equality.inertia={}/{}/{}",
        evidence.equality.inertia.positive,
        evidence.equality.inertia.negative,
        evidence.equality.inertia.zero
    );
    println!(
        "equality.backward_error={:.17e}",
        evidence.equality.normalized_backward_error
    );
    println!(
        "equality.recovery_round_trip={:.17e}",
        evidence.equality.scaling_round_trip_error
    );
    print_convex("qp", &evidence.qp);
    print_convex("socp", &evidence.socp);
    println!(
        "canonical.cross_route_error={:.17e}",
        evidence.cross_route_observable_error
    );
    println!(
        "canonical.field_energy={:.17e}",
        evidence.equality.recovered.field_energy
    );
    println!(
        "canonical.objective={:.17e}",
        evidence.equality.recovered.objective
    );
    println!(
        "counterexample.rank_deficient={:?}",
        counterexamples.rank_deficient_polynomial.kind
    );
    println!(
        "counterexample.nonpositive={:?}",
        counterexamples.nonpositive_reduced_pairing.kind
    );
    println!(
        "counterexample.broken_recovery={:?}",
        counterexamples.broken_recovery.kind
    );
    println!("counterexample.hidden_regularization=false");
    Ok(())
}

fn print_convex(label: &str, evidence: &georbf_cubic_cpd_recovery_spike::ConvexRouteEvidence) {
    println!("{label}.scaled.primal={:.17e}", evidence.scaled.primal);
    println!("{label}.scaled.dual={:.17e}", evidence.scaled.dual);
    println!(
        "{label}.scaled.stationarity={:.17e}",
        evidence.scaled.stationarity
    );
    println!(
        "{label}.scaled.complementarity={:.17e}",
        evidence.scaled.complementarity
    );
    println!(
        "{label}.scaled.relative_gap={:.17e}",
        evidence.scaled.relative_gap
    );
    println!(
        "{label}.recovery_round_trip={:.17e}",
        evidence.recovery_round_trip_error
    );
    println!(
        "{label}.physical.slack_equation_violation={:.17e}",
        evidence.physical_slack_equation_violation
    );
    println!(
        "{label}.canonical.hard_violation={:.17e}",
        evidence.recovered.hard_violation
    );
    println!(
        "{label}.canonical.truth_error={:.17e}",
        evidence.manufactured_truth_error
    );
}
