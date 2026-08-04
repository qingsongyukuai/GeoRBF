#[path = "../examples/convex_relations.rs"]
mod convex_relations_example;

#[test]
fn user_can_run_the_complete_v0_2_convex_relations_preview() {
    convex_relations_example::run().expect("the public v0.2.0 release example must run end to end");
}

#[test]
#[ignore = "the release workflow runs this quantity smoke with optimized code"]
fn v0_2_qp_smoke_accepts_512_constraints_and_10_000_ordered_queries() {
    convex_relations_example::run_smoke()
        .expect("the v0.2.0 QP quantity smoke must stay inside the accepted envelope");
}
