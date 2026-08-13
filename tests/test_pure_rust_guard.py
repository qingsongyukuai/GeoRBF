"""Tests for the repository's pure-Rust dependency and package guard."""

from __future__ import annotations

import importlib.util
import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
GUARD_PATH = ROOT / "tools" / "audit_pure_rust.py"


def load_guard():
    spec = importlib.util.spec_from_file_location("audit_pure_rust", GUARD_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load guard from {GUARD_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PureRustGuardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.guard = load_guard()

    def test_rejects_every_frozen_dependency_family(self) -> None:
        forbidden = {
            "cc": "native compiler helper",
            "cmake": "CMake",
            "bindgen": "bindgen",
            "cxx": "CXX bridge",
            "cxx-build": "CXX bridge",
            "eigen-sys": "Eigen",
            "qt_core": "Qt",
            "vtk-sys": "VTK",
            "pybind11-sys": "pybind11",
            "blas-src": "BLAS",
            "openblas-src": "OpenBLAS",
            "lapack-sys": "LAPACK",
            "intel-mkl-src": "MKL",
        }
        for package, category in forbidden.items():
            with self.subTest(package=package):
                self.assertEqual(
                    self.guard.forbidden_dependency_reason(package), category
                )

    def test_does_not_use_substring_matches_for_unrelated_crates(self) -> None:
        for package in ("success", "accounting", "black-box", "quote", "vtkio"):
            with self.subTest(package=package):
                self.assertIsNone(self.guard.forbidden_dependency_reason(package))

    def test_rejects_custom_build_targets_and_cargo_links(self) -> None:
        metadata = {
            "packages": [
                {
                    "name": "safe-rust",
                    "version": "1.0.0",
                    "links": None,
                    "targets": [{"kind": ["lib"], "name": "safe_rust"}],
                },
                {
                    "name": "generated-code",
                    "version": "1.0.0",
                    "links": None,
                    "targets": [{"kind": ["custom-build"], "name": "build-script-build"}],
                },
                {
                    "name": "native-wrapper",
                    "version": "1.0.0",
                    "links": "native_wrapper",
                    "targets": [{"kind": ["lib"], "name": "native_wrapper"}],
                },
            ]
        }
        violations = self.guard.metadata_violations(metadata)
        self.assertTrue(any("custom build target" in item for item in violations))
        self.assertTrue(any("Cargo links" in item for item in violations))

    def test_rejects_native_and_reference_paths(self) -> None:
        rejected = [
            "build.rs",
            "native/kernel.cpp",
            "native/kernel.hpp",
            "lib/georbf.so",
            "cmake/toolchain.cmake",
            ".cache/surfe-reference/CMakeLists.txt",
            ".cache/surfe-oracle/oracle",
        ]
        violations = self.guard.path_violations(rejected, "synthetic")
        self.assertEqual(len(violations), len(rejected))
        self.assertEqual(
            self.guard.path_violations(
                ["src/lib.rs", "docs/native-dependencies.md", "tests/data.json"],
                "synthetic",
            ),
            [],
        )

    def test_rejects_native_ffi_in_rust_sources(self) -> None:
        self.assertEqual(
            len(
                self.guard.rust_source_violations(
                    "src/native.rs",
                    '#[link(name = "blas")]\nunsafe extern "C" { fn dgemm(); }\n',
                    "synthetic",
                )
            ),
            2,
        )
        self.assertEqual(
            self.guard.rust_source_violations(
                "src/lib.rs", "pub fn rust_only() {}\n", "synthetic"
            ),
            [],
        )

    def test_current_repository_passes_the_guard(self) -> None:
        self.assertEqual(self.guard.audit_repository(ROOT), [])


if __name__ == "__main__":
    unittest.main()
