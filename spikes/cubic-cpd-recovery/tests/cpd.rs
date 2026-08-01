use georbf_cubic_cpd_recovery_spike::run_manufactured_experiment;

#[test]
fn manufactured_functionals_keep_the_complete_affine_polynomial() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");

    assert_eq!(evidence.cpd.polynomial_dimension, 4);
    assert_eq!(evidence.cpd.polynomial_rank, 4);
}

#[test]
fn cubic_pairing_is_strictly_positive_on_the_polynomial_null_space() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");

    assert!(evidence.cpd.null_space_defect <= evidence.acceptance.null_space_defect);
    assert!(evidence.cpd.reduced_symmetry_defect <= evidence.cpd.symmetry_defect_limit);
    assert!(evidence.cpd.reduced_smallest_eigenvalue > 0.0);
}

#[test]
fn complete_affine_fields_are_reproduced_without_kernel_energy() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");

    assert!(evidence.cpd.affine_reproduction_error <= evidence.acceptance.affine_reproduction);
}
