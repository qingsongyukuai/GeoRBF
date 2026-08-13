#!/usr/bin/env python3
"""Reject native dependencies, build scripts, and native package contents.

The guard uses Cargo's full all-feature metadata graph, so direct, transitive,
normal, development, build, optional, and target-specific packages are all
subject to the same policy. It uses only the Python standard library and never
discovers or starts the external Surfe reference or oracle.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys
from collections.abc import Iterable
from typing import Any


NATIVE_SUFFIXES = frozenset(
    {
        ".a",
        ".c",
        ".cc",
        ".cmake",
        ".cpp",
        ".cxx",
        ".dll",
        ".dylib",
        ".exe",
        ".h",
        ".hh",
        ".hpp",
        ".hxx",
        ".lib",
        ".o",
        ".obj",
        ".pdb",
        ".rlib",
        ".rmeta",
        ".so",
    }
)
REFERENCE_SEGMENTS = frozenset({"surfe-reference", "surfe-oracle"})
NATIVE_ABI_PATTERN = re.compile(
    r'\bextern\s*"(?:C|C-unwind|cdecl|stdcall|system|sysv64|win64)"'
)
LINK_ATTRIBUTE_PATTERN = re.compile(r"#\s*\[\s*link\s*\(")


def forbidden_dependency_reason(package_name: str) -> str | None:
    """Return the frozen-policy category rejected for a Cargo package name."""

    normalized = package_name.lower().replace("_", "-")
    tokens = tuple(token for token in normalized.split("-") if token)

    if normalized in {"cc", "gcc"}:
        return "native compiler helper"
    if "cmake" in tokens:
        return "CMake"
    if any(token.startswith("bindgen") for token in tokens):
        return "bindgen"
    if normalized == "cxx" or normalized.startswith(("cxx-", "cxxbridge-")):
        return "CXX bridge"
    if "eigen" in tokens:
        return "Eigen"
    if any(re.fullmatch(r"qt[0-9]*", token) for token in tokens) or normalized in {
        "qmetaobject",
        "ritual",
    }:
        return "Qt"
    if normalized == "vtk" or normalized.startswith("vtk-"):
        return "VTK"
    if any(token.startswith("pybind") for token in tokens):
        return "pybind11"
    if "openblas" in normalized:
        return "OpenBLAS"
    if "lapack" in tokens or normalized.startswith("lapack-"):
        return "LAPACK"
    if "mkl" in tokens or "intel-mkl" in normalized:
        return "MKL"
    if "blas" in tokens or normalized.startswith("blas-"):
        return "BLAS"
    return None


def metadata_violations(metadata: dict[str, Any]) -> list[str]:
    """Audit every package represented in `cargo metadata` output."""

    violations: list[str] = []
    for package in metadata.get("packages", []):
        name = str(package.get("name", "<unnamed>"))
        version = str(package.get("version", "<unknown>"))
        label = f"{name}@{version}"
        if reason := forbidden_dependency_reason(name):
            violations.append(f"dependency {label}: forbidden {reason} package")
        if links := package.get("links"):
            violations.append(f"dependency {label}: Cargo links={links!r} is forbidden")
        for target in package.get("targets", []):
            if "custom-build" in target.get("kind", []):
                violations.append(f"dependency {label}: custom build target is forbidden")
                break
    return sorted(set(violations))


def path_violations(paths: Iterable[str], source: str) -> list[str]:
    """Audit tracked or packaged paths for native and oracle content."""

    violations: list[str] = []
    for raw_path in paths:
        path = raw_path.strip().replace("\\", "/").removeprefix("./")
        if not path:
            continue
        pure_path = pathlib.PurePosixPath(path)
        lower_parts = tuple(part.lower() for part in pure_path.parts)
        lower_name = pure_path.name.lower()
        suffix = pure_path.suffix.lower()
        reason: str | None = None
        if lower_name == "build.rs":
            reason = "build script"
        elif lower_name == "cmakelists.txt" or suffix == ".cmake":
            reason = "CMake input"
        elif suffix in NATIVE_SUFFIXES:
            reason = "native source or compiled artifact"
        elif REFERENCE_SEGMENTS.intersection(lower_parts):
            reason = "Surfe reference/oracle content"
        if reason:
            violations.append(f"{source} path {path!r}: forbidden {reason}")
    return violations


def rust_source_violations(path: str, contents: str, source: str) -> list[str]:
    """Reject explicit native ABI and linker declarations in Rust source."""

    violations: list[str] = []
    if NATIVE_ABI_PATTERN.search(contents):
        violations.append(f"{source} Rust source {path!r}: native extern ABI is forbidden")
    if LINK_ATTRIBUTE_PATTERN.search(contents):
        violations.append(f"{source} Rust source {path!r}: #[link] is forbidden")
    return violations


def metadata_rust_source_violations(metadata: dict[str, Any]) -> list[str]:
    """Scan Rust sources for every direct and transitive Cargo package."""

    violations: list[str] = []
    scanned: set[pathlib.Path] = set()
    for package in metadata.get("packages", []):
        manifest = pathlib.Path(str(package.get("manifest_path", "")))
        if not manifest.is_file():
            continue
        package_root = manifest.parent
        label = f"dependency {package.get('name', '<unnamed>')}@{package.get('version', '<unknown>')}"
        for source_path in package_root.rglob("*.rs"):
            if source_path in scanned or any(
                part in {".git", ".cache", "target"} for part in source_path.parts
            ):
                continue
            scanned.add(source_path)
            try:
                contents = source_path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                violations.append(f"{label} Rust source {source_path}: invalid UTF-8")
                continue
            violations.extend(rust_source_violations(str(source_path), contents, label))
    return violations


def run(command: list[str], root: pathlib.Path) -> tuple[str, str | None]:
    """Run a read-only audit command and return stdout or a stable error."""

    completed = subprocess.run(
        command,
        cwd=root,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        return "", f"command {' '.join(command)!r} failed: {detail}"
    return completed.stdout, None


def manifest_violations(root: pathlib.Path) -> list[str]:
    """Reject explicit Cargo `links` declarations even outside the resolve graph."""

    violations: list[str] = []
    for manifest in root.rglob("Cargo.toml"):
        relative = manifest.relative_to(root)
        if any(part in {".git", ".cache", "target"} for part in relative.parts):
            continue
        contents = manifest.read_text(encoding="utf-8")
        if re.search(r"(?m)^\s*links\s*=", contents):
            violations.append(f"manifest {relative}: Cargo links declaration is forbidden")
    return violations


def audit_repository(root: pathlib.Path) -> list[str]:
    """Run the complete local/CI pure-Rust audit for a repository."""

    root = root.resolve()
    violations: list[str] = []

    metadata_output, error = run(
        [
            "cargo",
            "metadata",
            "--locked",
            "--all-features",
            "--format-version",
            "1",
        ],
        root,
    )
    if error:
        violations.append(error)
    else:
        try:
            metadata = json.loads(metadata_output)
            violations.extend(metadata_violations(metadata))
            violations.extend(metadata_rust_source_violations(metadata))
        except json.JSONDecodeError as exception:
            violations.append(f"cargo metadata returned invalid JSON: {exception}")

    repository_paths, error = run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard"], root
    )
    if error:
        violations.append(error)
    else:
        violations.extend(path_violations(repository_paths.splitlines(), "repository"))

    package_paths, error = run(
        ["cargo", "package", "--list", "--locked", "--allow-dirty"], root
    )
    if error:
        violations.append(error)
    else:
        violations.extend(path_violations(package_paths.splitlines(), "package"))

    violations.extend(manifest_violations(root))
    return sorted(set(violations))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parents[1],
        help="repository root (defaults to the guard's parent repository)",
    )
    arguments = parser.parse_args(argv)
    violations = audit_repository(arguments.root)
    if violations:
        print("pure-Rust audit failed:", file=sys.stderr)
        for violation in violations:
            print(f"- {violation}", file=sys.stderr)
        return 1
    print("pure-Rust audit passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
