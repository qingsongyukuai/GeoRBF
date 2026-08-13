use georbf::{
    AnisotropyError, BuildError, Builder, Constraints, ContinuousPropertyError, DenseMatrix, Error,
    Inequality, Interface, ModelType, Planar, RbfKernel,
};

#[test]
fn stable_surfe_error_categories_remain_exhaustive() {
    assert_eq!(Error::ALL.len(), 23);
    for category in Error::ALL {
        assert!(!category.surfe_exception_name().is_empty());
        assert!(!category.message().is_empty());
        assert_eq!(category.to_string(), category.message());
    }
}

#[test]
fn public_rejection_matrix_uses_categories_not_message_matching() {
    assert_eq!(
        Builder::from_surfe_model_code(0).unwrap_err(),
        BuildError::Surfe(Error::UnknownModel)
    );

    let mut unknown_kernel = Builder::new(ModelType::SingleSurface);
    assert_eq!(
        unknown_kernel.set_rbf_kernel_name("cubic").unwrap_err(),
        BuildError::Surfe(Error::UnknownRbf)
    );

    let mut no_interface = Builder::new(ModelType::SingleSurface);
    no_interface
        .add_planar_normal(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
        .unwrap();
    assert_eq!(
        no_interface.fit().unwrap_err().surfe_category(),
        Some(Error::NoInterfaceData)
    );

    let one_interface = Builder::new(ModelType::LajaunieApproach)
        .set_constraints(georbf::Constraints {
            interfaces: vec![Interface::new(0.0, 0.0, 0.0, 1.0).unwrap()],
            ..georbf::Constraints::default()
        })
        .fit()
        .unwrap_err();
    assert_eq!(
        one_interface.surfe_category(),
        Some(Error::NoInterfaceIncrementPairs)
    );

    let invalid_stratigraphy = Builder::new(ModelType::StratigraphicHorizons)
        .set_constraints(Constraints {
            interfaces: vec![
                Interface::new(0.0, 0.0, 0.0, 1.0).unwrap(),
                Interface::new(1.0, 0.0, 0.0, 1.0).unwrap(),
            ],
            inequalities: vec![Inequality::new(0.0, 1.0, 0.0, 1.0).unwrap()],
            ..Constraints::default()
        })
        .fit()
        .unwrap_err();
    assert_eq!(
        invalid_stratigraphy.surfe_category(),
        Some(Error::InvalidInputData)
    );

    let wrong_shape = DenseMatrix::from_row_major(1, 3, vec![0.0; 3]).unwrap();
    assert_eq!(
        Builder::new(ModelType::SingleSurface)
            .set_interface_constraint_matrix(&wrong_shape)
            .unwrap_err()
            .surfe_category(),
        Some(Error::IncorrectArrayDimensions)
    );

    let fitted = Builder::new(ModelType::SingleSurface)
        .set_constraints(Constraints {
            interfaces: vec![
                Interface::new(0.0, 0.0, 0.0, 0.0).unwrap(),
                Interface::new(1.0, 0.0, 0.0, 0.0).unwrap(),
                Interface::new(0.0, 1.0, 0.0, 0.0).unwrap(),
                Interface::new(1.0, 1.0, 0.0, 0.0).unwrap(),
            ],
            planars: vec![Planar::from_normal(0.5, 0.5, 0.0, 0.0, 0.0, 1.0).unwrap()],
            ..Constraints::default()
        })
        .fit()
        .unwrap();
    let wrong_evaluation_shape = DenseMatrix::from_row_major(1, 2, vec![0.0; 2]).unwrap();
    assert_eq!(
        fitted
            .evaluate_scalar_matrix(&wrong_evaluation_shape)
            .unwrap_err()
            .surfe_category(),
        Some(Error::IncorrectArrayDimensions)
    );
}

#[test]
fn setup_and_solver_failures_map_to_the_frozen_outer_category() {
    let mut anisotropy = Builder::new(ModelType::ContinuousProperty);
    anisotropy
        .add_interface(Interface::new(0.0, 0.0, 0.0, 1.0).unwrap())
        .set_global_anisotropy(true);
    assert_eq!(
        anisotropy.fit().unwrap_err().surfe_category(),
        Some(Error::BasisFunctionSetupFailure)
    );

    let mut unsupported_anisotropy = Builder::new(ModelType::ContinuousProperty);
    unsupported_anisotropy
        .add_interface(Interface::new(0.0, 0.0, 0.0, 1.0).unwrap())
        .add_interface(Interface::new(1.0, 0.0, 0.0, 2.0).unwrap())
        .add_planar(Planar::from_normal(0.0, 0.0, 0.0, 1.0, 0.0, 0.0).unwrap())
        .add_planar(Planar::from_normal(1.0, 1.0, 1.0, 0.0, 1.0, 0.0).unwrap())
        .set_rbf_kernel(RbfKernel::WendlandC2)
        .set_global_anisotropy(true);
    assert_eq!(
        unsupported_anisotropy.fit().unwrap_err().surfe_category(),
        Some(Error::BasisFunctionSetupFailure)
    );

    let mut singular = Builder::new(ModelType::ContinuousProperty);
    singular.add_interface(Interface::new(0.0, 0.0, 0.0, 1.0).unwrap());
    assert_eq!(
        singular.fit().unwrap_err().surfe_category(),
        Some(Error::LinearSolverFailure)
    );
}

#[test]
fn safe_rust_rejections_are_documented_as_no_direct_surfe_category() {
    let mut non_finite = Builder::new(ModelType::SingleSurface);
    non_finite.set_rbf_shape_parameter(f64::NAN);
    assert_eq!(non_finite.fit().unwrap_err().surfe_category(), None);

    let safe = BuildError::ContinuousProperty(ContinuousPropertyError::EqualityVectorOutOfBounds {
        planar_count: 1,
        tangent_count: 0,
    });
    assert_eq!(safe.surfe_category(), None);

    let nested = BuildError::ContinuousProperty(ContinuousPropertyError::Anisotropy(
        AnisotropyError::InsufficientPlanars,
    ));
    assert_eq!(
        nested.surfe_category(),
        Some(Error::BasisFunctionSetupFailure)
    );
}
