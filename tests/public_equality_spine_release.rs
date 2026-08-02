#[path = "../examples/equality_spine.rs"]
mod equality_spine_example;

#[test]
fn user_can_run_the_complete_v0_1_equality_spine() {
    equality_spine_example::run().expect("the public release example must run end to end");
}
