"""Validate release artifacts that are independent of requirement tracing."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path


PLACEHOLDER_PATTERNS = (
    re.compile(r"\btodo!\s*\("),
    re.compile(r"\bunimplemented!\s*\("),
    re.compile(r"not[ _-]implemented", re.IGNORECASE),
)


def audit_source_placeholders(root: Path) -> list[str]:
    failures: list[str] = []
    source_root = root / "src"
    if not source_root.is_dir():
        return ["missing product source directory src"]
    for path in sorted(source_root.rglob("*.rs")):
        content = path.read_text(encoding="utf-8")
        for pattern in PLACEHOLDER_PATTERNS:
            match = pattern.search(content)
            if match is not None:
                line = content.count("\n", 0, match.start()) + 1
                failures.append(
                    f"unsupported implementation placeholder in "
                    f"{path.relative_to(root).as_posix()}:{line}: {match.group(0)}"
                )
    return failures


def mirrored_files(root: Path, subtree: str) -> set[Path]:
    directory = root / subtree
    if not directory.is_dir():
        return set()
    return {
        path.relative_to(directory)
        for path in directory.rglob("*")
        if path.is_file()
    }


def audit_oracle_mirror(root: Path) -> list[str]:
    failures: list[str] = []
    spike = root / "spikes/oracle-fixtures"
    adopted = root / "validation/oracle/cubic-v1"
    pairs = [
        (spike / "manifest.json", adopted / "source-manifest.json", Path("manifest.json"))
    ]
    for subtree in ("cases", "fixtures"):
        spike_files = mirrored_files(spike, subtree)
        adopted_files = mirrored_files(adopted, subtree)
        if spike_files != adopted_files:
            missing = sorted(path.as_posix() for path in spike_files - adopted_files)
            unexpected = sorted(path.as_posix() for path in adopted_files - spike_files)
            if missing:
                failures.append(f"oracle adoption missing {subtree}: {', '.join(missing)}")
            if unexpected:
                failures.append(
                    f"oracle adoption has unexpected {subtree}: {', '.join(unexpected)}"
                )
        for relative in sorted(spike_files & adopted_files):
            pairs.append(
                (spike / subtree / relative, adopted / subtree / relative, Path(subtree) / relative)
            )
    for source, destination, display in pairs:
        if not source.is_file() or not destination.is_file():
            failures.append(f"oracle mirror path is missing: {display.as_posix()}")
        elif source.read_bytes() != destination.read_bytes():
            failures.append(f"oracle adoption drift: {display.as_posix()}")
    return failures


def require_markers(root: Path, path_value: str, markers: tuple[str, ...]) -> list[str]:
    path = root / path_value
    if not path.is_file():
        return [f"missing release artifact: {path_value}"]
    try:
        content = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        return [f"cannot read release artifact {path_value}: {error}"]
    return [
        f"release artifact {path_value} is missing marker {marker!r}"
        for marker in markers
        if marker not in content
    ]


def audit_release_metadata(root: Path) -> list[str]:
    failures: list[str] = []
    manifest_path = root / "Cargo.toml"
    if not manifest_path.is_file():
        return ["missing Cargo.toml"]
    try:
        manifest = tomllib.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as error:
        return [f"cannot parse Cargo.toml: {error}"]
    package = manifest.get("package", {})
    expected_fields = {"version": "0.2.0", "edition": "2024", "rust-version": "1.85"}
    for field, expected in expected_fields.items():
        if package.get(field) != expected:
            failures.append(f"Cargo package.{field} must be {expected}")
    expected_package_paths = {
        "/RELEASE_NOTES.md",
        "/examples/**",
        "/docs/implementation/25-equality-spine-release.md",
        "/docs/implementation/36-convex-relations-release.md",
        "/validation/v0.1.0/**",
        "/validation/v0.2.0/**",
    }
    included = set(package.get("include", []))
    missing_package_paths = sorted(expected_package_paths - included)
    if missing_package_paths:
        failures.append(
            "Cargo package excludes release artifacts: " + ", ".join(missing_package_paths)
        )

    artifact_markers = {
        "RELEASE_NOTES.md": (
            "# GeoRBF v0.2.0",
            "## Supported scope",
            "## Compatibility boundary",
            "## Diagnostic semantics",
            "## Out of scope",
            "## Verification",
        ),
        "docs/implementation/36-convex-relations-release.md": (
            "Issue: [#36]",
            "Evidence seams: T01–T17",
            "Traceability audit",
            "Release procedure",
        ),
        "examples/convex_relations.rs": (
            "pub fn run",
            "HorizonBuilder",
            "CovarianceGroupBuilder",
            "DirectionalDerivativeInterval",
            "PolarityResolution",
            "evaluate_batch",
            "run_smoke",
        ),
        "README.md": ("Version 0.2.0", "examples/convex_relations.rs"),
        ".github/workflows/product-v0.2.yml": (
            "v0.2.0",
            "release-corpus",
            "PROPTEST_CASES",
            "PROPTEST_RNG_SEED",
            "scripts/release_audit.py",
            "spikes/oracle-fixtures/generate.py --check",
            "cargo build --locked --release",
            "v0_2_qp_smoke",
            "--example convex_relations",
            "cargo package --locked",
        ),
    }
    for path, markers in artifact_markers.items():
        failures.extend(require_markers(root, path, markers))
    return failures


def audit_repository_artifacts(root: Path) -> list[str]:
    failures = audit_source_placeholders(root)
    failures.extend(audit_oracle_mirror(root))
    failures.extend(audit_release_metadata(root))
    return failures
