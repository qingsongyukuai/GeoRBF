import tempfile
import unittest
from pathlib import Path

from release_contract import RELEASE_VERSION
from release_repository_checks import audit_oracle_mirror, audit_source_placeholders
from release_traceability import EXPECTED_REQUIREMENTS, validate_traceability


class TraceabilityAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        (self.root / "tests").mkdir()
        (self.root / "docs").mkdir()
        (self.root / "tests" / "evidence.rs").write_text(
            "fn cumulative_release_evidence() {}\n", encoding="utf-8"
        )
        (self.root / "docs" / "release.md").write_text(
            "# Release evidence\n", encoding="utf-8"
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def traceability(self) -> dict[str, object]:
        return {
            "schema_version": "georbf-traceability-v1",
            "release": RELEASE_VERSION,
            "evidence_sets": {
                "release": [
                    {
                        "path": "tests/evidence.rs",
                        "contains": "cumulative_release_evidence",
                    }
                ]
            },
            "documentation_sets": {
                "release": [
                    {"path": "docs/release.md", "contains": "Release evidence"}
                ]
            },
            "requirements": [
                {
                    "id": requirement,
                    "behavior": f"Release behavior for {requirement}.",
                    "api_ir": f"Public API or IR path for {requirement}",
                    "evidence": "release",
                    "documentation": "release",
                }
                for requirement in sorted(EXPECTED_REQUIREMENTS)
            ],
        }

    def test_complete_traceability_is_accepted(self) -> None:
        self.assertEqual(validate_traceability(self.root, self.traceability()), [])

    def test_v02_requirement_scope_matches_the_convex_relations_milestone(self) -> None:
        self.assertIn("PAPI-019", EXPECTED_REQUIREMENTS)
        self.assertIn("DOM-022", EXPECTED_REQUIREMENTS)
        self.assertIn("VAL-015", EXPECTED_REQUIREMENTS)
        self.assertNotIn("PAPI-016", EXPECTED_REQUIREMENTS)
        self.assertNotIn("DOM-008", EXPECTED_REQUIREMENTS)
        self.assertNotIn("KER-002", EXPECTED_REQUIREMENTS)
        self.assertNotIn("VAL-016", EXPECTED_REQUIREMENTS)

    def test_duplicate_and_missing_requirements_fail_closed(self) -> None:
        traceability = self.traceability()
        requirements = traceability["requirements"]
        requirements[1]["id"] = requirements[0]["id"]

        failures = validate_traceability(self.root, traceability)

        self.assertTrue(any("duplicate requirement" in failure for failure in failures))
        self.assertTrue(any("missing requirements" in failure for failure in failures))

    def test_dangling_evidence_and_documentation_fail_closed(self) -> None:
        traceability = self.traceability()
        traceability["evidence_sets"]["release"][0]["contains"] = "missing-test"
        traceability["documentation_sets"]["release"][0]["path"] = (
            "docs/missing.md"
        )

        failures = validate_traceability(self.root, traceability)

        self.assertTrue(any("missing marker" in failure for failure in failures))
        self.assertTrue(any("missing referenced path" in failure for failure in failures))

    def test_duplicate_behavior_and_api_ir_paths_fail_closed(self) -> None:
        traceability = self.traceability()
        requirements = traceability["requirements"]
        requirements[1]["behavior"] = requirements[0]["behavior"]
        requirements[1]["api_ir"] = requirements[0]["api_ir"]

        failures = validate_traceability(self.root, traceability)

        self.assertTrue(any("behaviors must be unique" in failure for failure in failures))
        self.assertTrue(any("API/IR paths must be unique" in failure for failure in failures))

    def test_unknown_and_unused_reference_sets_fail_closed(self) -> None:
        traceability = self.traceability()
        traceability["evidence_sets"]["unused"] = [
            {
                "path": "tests/evidence.rs",
                "contains": "cumulative_release_evidence",
            }
        ]
        traceability["requirements"][0]["documentation"] = "missing"

        failures = validate_traceability(self.root, traceability)

        self.assertTrue(any("names unknown set" in failure for failure in failures))
        self.assertTrue(any("unused evidence_sets" in failure for failure in failures))


class ReleaseSourceAuditTests(unittest.TestCase):
    def test_product_placeholders_are_rejected_without_scanning_tests(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "src").mkdir()
            (root / "tests").mkdir()
            (root / "src" / "lib.rs").write_text(
                'fn unsupported() { todo!("later"); }\n', encoding="utf-8"
            )
            (root / "tests" / "fixture.rs").write_text(
                'const MESSAGE: &str = "not implemented";\n', encoding="utf-8"
            )

            failures = audit_source_placeholders(root)

            self.assertEqual(len(failures), 1)
            self.assertIn("src/lib.rs", failures[0])

    def test_oracle_adoption_must_be_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            spike = root / "spikes" / "oracle-fixtures"
            adopted = root / "validation" / "oracle" / "cubic-v1"
            for relative in ("cases", "fixtures"):
                (spike / relative).mkdir(parents=True)
                (adopted / relative).mkdir(parents=True)
                (spike / relative / "case.json").write_text("{}\n", encoding="utf-8")
                (adopted / relative / "case.json").write_text("{}\n", encoding="utf-8")
            (spike / "manifest.json").write_text("{}\n", encoding="utf-8")
            (adopted / "source-manifest.json").write_text("{}\n", encoding="utf-8")

            self.assertEqual(audit_oracle_mirror(root), [])

            (adopted / "fixtures" / "case.json").write_text(
                '{"drift": true}\n', encoding="utf-8"
            )
            failures = audit_oracle_mirror(root)
            self.assertEqual(len(failures), 1)
            self.assertIn("fixtures/case.json", failures[0])


if __name__ == "__main__":
    unittest.main()
