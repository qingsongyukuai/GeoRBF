from __future__ import annotations

import pathlib
import tomllib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "290dbe0ab344f4258a4935f05cad0f153f0f69a4"


class ReleaseAuditTests(unittest.TestCase):
    def test_cargo_release_metadata_and_package_documents(self) -> None:
        manifest = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
        package = manifest["package"]

        self.assertEqual(package["license"], "MIT")
        self.assertEqual(package["readme"], "README.md")
        included = set(package["include"])
        for required in {
            "/LICENSE",
            "/NOTICE",
            "/README.md",
            "/docs/port/compatibility.md",
            "/docs/port/licensing-and-rust-boundary.md",
            "/docs/port/parity-report.md",
            "/docs/port/performance-report.md",
            "/docs/port/release-audit.md",
            "/docs/port/source-traceability.md",
        }:
            self.assertIn(required, included)

    def test_upstream_notice_is_complete_and_pinned(self) -> None:
        notice = (ROOT / "NOTICE").read_text(encoding="utf-8")
        for required in {
            "Copyright (c) 2017 Government of Canada",
            SOURCE_COMMIT,
            "MIT License",
            "Permission is hereby granted, free of charge",
            'THE SOFTWARE IS PROVIDED "AS IS"',
        }:
            self.assertIn(required, notice)

    def test_readme_documents_safe_public_lifecycle_and_evidence(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        for required in {
            "pure Rust",
            "Builder",
            "FittedModel",
            "SingleSurface",
            "LajaunieApproach",
            "StratigraphicHorizons",
            "ContinuousProperty",
            "VectorField",
            "docs/port/compatibility.md",
            "docs/port/parity-report.md",
            "docs/port/performance-report.md",
            "docs/port/release-audit.md",
            "docs/port/source-traceability.md",
        }:
            self.assertIn(required, readme)

    def test_every_production_module_has_pinned_source_traceability(self) -> None:
        traceability = (ROOT / "docs/port/source-traceability.md").read_text(
            encoding="utf-8"
        )
        lines = traceability.splitlines()
        production = sorted(
            path.relative_to(ROOT).as_posix() for path in (ROOT / "src").rglob("*.rs")
        )

        self.assertGreater(len(production), 30)
        for module in production:
            evidence = [line for line in lines if f"`{module}`" in line]
            self.assertEqual(len(evidence), 1, module)
            self.assertIn(f"@{SOURCE_COMMIT}", evidence[0], module)

    def test_final_audit_links_all_four_release_gates(self) -> None:
        audit = (ROOT / "docs/port/release-audit.md").read_text(encoding="utf-8")
        for required in {
            SOURCE_COMMIT,
            "behavior_parity",
            "pure_rust",
            "performance_not_lower_than_surfe",
            "release_audit",
            "parity-report.md",
            "performance-report.md",
            "source-traceability.md",
            "Linux",
            "macOS",
            "Windows",
        }:
            self.assertIn(required, audit)

    def test_ci_matrix_runs_the_release_blocking_commands(self) -> None:
        workflow = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")
        for required in {
            "workflow_dispatch:",
            "ubuntu-latest",
            "macos-latest",
            "windows-latest",
            "cargo fmt --all --check",
            "cargo clippy --all-targets --all-features -- -D warnings",
            "cargo test --all-targets --all-features",
            "cargo test --doc --all-features",
            "cargo build --release --all-features",
            "cargo package --locked",
            "python tools/audit_pure_rust.py",
            "tests/test_release_audit.py",
        }:
            self.assertIn(required, workflow)


if __name__ == "__main__":
    unittest.main()
