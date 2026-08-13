//! Release-blocking crosswalk for the formal T32 Surfe parity fixture.

#[path = "../common/parity/mod.rs"]
pub mod fixture;

use fixture::{parse_json, sha256_hex, validate_fixture, JsonValue};
use std::collections::BTreeSet;
use std::path::Path;

const GLOBAL_FIXTURE: &str = include_str!("../fixtures/golden/global-parity-v1.json");

fn object(value: &JsonValue) -> &[(String, JsonValue)] {
    match value {
        JsonValue::Object(fields) => fields,
        other => panic!("expected object, got {other:?}"),
    }
}

fn field<'a>(value: &'a JsonValue, name: &str) -> &'a JsonValue {
    object(value)
        .iter()
        .find_map(|(candidate, value)| (candidate == name).then_some(value))
        .unwrap_or_else(|| panic!("missing field {name}"))
}

fn string<'a>(value: &'a JsonValue, name: &str) -> &'a str {
    match field(value, name) {
        JsonValue::String(value) => value,
        other => panic!("expected string field {name}, got {other:?}"),
    }
}

fn number(value: &JsonValue, name: &str) -> usize {
    match field(value, name) {
        JsonValue::Number(value) => value
            .parse()
            .unwrap_or_else(|_| panic!("field {name} is not an integer")),
        other => panic!("expected number field {name}, got {other:?}"),
    }
}

fn array<'a>(value: &'a JsonValue, name: &str) -> &'a [JsonValue] {
    match field(value, name) {
        JsonValue::Array(values) => values,
        other => panic!("expected array field {name}, got {other:?}"),
    }
}

fn is_sha256(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn parsed_fixture() -> JsonValue {
    parse_json(GLOBAL_FIXTURE.trim_end()).expect("formal global fixture must be valid JSON")
}

fn capabilities(fixture: &JsonValue) -> &JsonValue {
    field(field(field(fixture, "expected"), "result"), "capabilities")
}

#[test]
fn formal_fixture_is_canonical_schema_valid_and_hash_bound() {
    let parsed = parsed_fixture();
    assert_eq!(validate_fixture(&parsed), Ok(()));
    assert_eq!(format!("{}\n", parsed.canonical_json()), GLOBAL_FIXTURE);

    let request = field(&parsed, "request");
    let expected = field(&parsed, "expected");
    let request_line = format!("{}\n", request.canonical_json());
    let response_line = format!("{}\n", expected.canonical_json());
    let dataset = field(&parsed, "dataset");
    assert_eq!(
        sha256_hex(request_line.as_bytes()),
        string(dataset, "request_line_sha256")
    );
    assert_eq!(
        sha256_hex(response_line.as_bytes()),
        string(dataset, "response_line_sha256")
    );
}

#[test]
fn every_frozen_probe_and_fixed_family_has_a_live_rust_parity_anchor() {
    let parsed = parsed_fixture();
    let capabilities = capabilities(&parsed);
    assert_eq!(string(capabilities, "suite"), "global-parity-v1");
    assert_eq!(
        field(capabilities, "reference_clean"),
        &JsonValue::Bool(true)
    );

    let probe_catalog = array(capabilities, "probe_catalog");
    let expected_tasks: BTreeSet<_> = (6..=31).map(|task| format!("T{task:02}")).collect();
    let mut actual_tasks = BTreeSet::new();
    for probe in probe_catalog {
        let task = string(probe, "task").to_owned();
        assert!(actual_tasks.insert(task), "duplicate probe task");
        for name in ["source_sha256", "binary_sha256", "transcript_sha256"] {
            assert!(is_sha256(string(probe, name)), "invalid {name}");
        }
        assert_eq!(
            sha256_hex(string(probe, "transcript").as_bytes()),
            string(probe, "transcript_sha256"),
            "transcript hash mismatch for {}",
            string(probe, "task")
        );
        assert!(!array(probe, "invocations").is_empty());
    }
    assert_eq!(actual_tasks, expected_tasks);

    let families = array(capabilities, "family_coverage");
    let mut family_ids = BTreeSet::new();
    let manifest_root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for row in families {
        let family = string(row, "family");
        assert!(
            family_ids.insert(family.to_owned()),
            "duplicate family {family}"
        );
        assert_eq!(string(row, "status"), "passed");
        assert!(actual_tasks.contains(string(row, "probe_task")));

        let anchor = string(row, "rust_test");
        let (relative_path, test_name) = anchor
            .rsplit_once("::")
            .unwrap_or_else(|| panic!("invalid Rust test anchor {anchor}"));
        let source = std::fs::read_to_string(manifest_root.join(relative_path))
            .unwrap_or_else(|error| panic!("cannot read {relative_path}: {error}"));
        assert!(
            source.contains(&format!("fn {test_name}(")),
            "missing Rust test anchor {anchor}"
        );
    }

    for kernel in [
        "cubic",
        "gaussian",
        "mq",
        "mq3",
        "tps",
        "imq",
        "r",
        "wendland_c2",
        "matern_c4",
    ] {
        for family in [
            format!("kernel/isotropic/{kernel}/separated"),
            format!("kernel/isotropic/{kernel}/zero-near-support"),
            format!("kernel/functionals/{kernel}/directions"),
            format!("kernel/modified/{kernel}/all-functionals"),
        ] {
            assert!(
                family_ids.contains(&family),
                "missing fixed family {family}"
            );
        }
    }
    for family in [
        "kernel/anisotropy/identity-oblique-degenerate",
        "model/single_surface/equality",
        "model/single_surface/inequality-active-inactive",
        "model/single_surface/restricted-range",
        "model/lajaunie/multilevel-increments",
        "model/stratigraphic/three-level-lithology",
        "model/continuous_property/reachable",
        "model/vector_field/planar-hessian",
        "solver/lu/well-conditioned",
        "solver/lu/ill-conditioned-attempted",
        "solver/lu/singular",
        "solver/lu/non-finite",
        "solver/qp/equality-inequality",
        "solver/qp/active-boundary",
        "solver/qp/inactive-boundary",
        "solver/qp/infeasible",
        "solver/loqo/single-double-tight-bound",
        "solver/loqo/infeasible",
    ] {
        assert!(family_ids.contains(family), "missing fixed family {family}");
    }
    for stage in [
        "request",
        "configuration",
        "constraint-ingest",
        "preprocess",
        "basis",
        "assembly",
        "solve",
        "reconstruction",
        "evaluation",
        "oracle-safety",
    ] {
        let family = format!("error/{stage}/typed-category");
        assert!(
            family_ids.contains(&family),
            "missing error family {family}"
        );
    }

    let summary = field(capabilities, "summary");
    assert_eq!(number(summary, "probe_tasks"), probe_catalog.len());
    assert_eq!(number(summary, "families"), families.len());
    assert_eq!(number(summary, "rust_tests_at_generation"), 198);
    assert_eq!(number(summary, "ignored_rust_tests"), 0);
    assert_eq!(number(summary, "unexplained_mismatches"), 0);
}
