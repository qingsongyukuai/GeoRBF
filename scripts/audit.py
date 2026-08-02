#!/usr/bin/env python3
"""Fail-closed audit of the product crate's active dependency envelope."""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


ALLOWED_LICENSE_EXPRESSIONS = {
    "(MIT OR Apache-2.0) AND Unicode-3.0",
    "Apache-2.0",
    "Apache-2.0 / MIT",
    "Apache-2.0 OR MIT",
    "BSD-2-Clause OR Apache-2.0 OR MIT",
    "BSD-3-Clause",
    "MIT",
    "MIT OR Apache-2.0",
    "MIT/Apache-2.0",
    "Unlicense OR MIT",
    "Unlicense/MIT",
    "Zlib OR Apache-2.0 OR MIT",
}
FORBIDDEN_PACKAGES = {
    "bindgen",
    "blas",
    "blas-src",
    "cc",
    "clang-sys",
    "cmake",
    "intel-mkl-src",
    "lapack",
    "lapack-src",
    "netlib-src",
    "openblas-src",
    "pardiso-wrapper",
    "pkg-config",
    "vcpkg",
}
FORBIDDEN_FEATURE_FRAGMENTS = (
    "blas",
    "lapack",
    "mkl",
    "netlib",
    "openblas",
    "pardiso",
    "sdp",
)
AUDITED_LOCKFILE_SHA256 = "1b6107e4b044251da66dea980740c0f65e38f1890dbd9211e3dd1ee1c938bb70"
AUDITED_BUILD_SCRIPTS = {
    "crunchy 0.2.4",
    "libc 0.2.189",
    "libm 0.2.16",
    "nano-gemm-c32 0.2.1",
    "nano-gemm-c64 0.2.1",
    "nano-gemm-f32 0.2.1",
    "nano-gemm-f64 0.2.1",
    "num-traits 0.2.19",
    "paste 1.0.15",
    "private-gemm-x86 0.1.20",
    "proc-macro2 1.0.107",
    "pulp 0.22.3",
    "quote 1.0.47",
    "syn 1.0.109",
    "thiserror 1.0.69",
    "zerocopy 0.8.55",
}


def run(*command: str) -> str:
    return subprocess.run(command, check=True, text=True, capture_output=True).stdout


def canonical_text_sha256(path: Path) -> str:
    """Hash UTF-8 text after normalizing platform-specific line endings."""
    text = path.read_text(encoding="utf-8")
    canonical_text = text.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(canonical_text.encode("utf-8")).hexdigest()


def main() -> int:
    rustc_verbose = run("rustc", "--version", "--verbose")
    rustc_version = rustc_verbose.splitlines()[0]
    host = next(
        line.removeprefix("host: ")
        for line in rustc_verbose.splitlines()
        if line.startswith("host: ")
    )
    expected_target = os.environ.get("EXPECTED_TARGET")
    if expected_target and expected_target != host:
        raise RuntimeError(
            f"runner target mismatch: expected {expected_target}, found {host}"
        )
    expected_rustc = os.environ.get("EXPECTED_RUSTC", "rustc 1.85.0")
    if rustc_version != expected_rustc and not rustc_version.startswith(
        f"{expected_rustc} "
    ):
        raise RuntimeError(
            f"toolchain mismatch: expected {expected_rustc}, found {rustc_version}"
        )

    metadata = json.loads(
        run(
            "cargo",
            "metadata",
            "--locked",
            "--format-version",
            "1",
            "--filter-platform",
            host,
        )
    )
    package_by_id = {package["id"]: package for package in metadata["packages"]}
    selected_nodes = metadata["resolve"]["nodes"]
    selected_packages = [package_by_id[node["id"]] for node in selected_nodes]
    selected_features = {
        package_by_id[node["id"]]["name"]: sorted(node["features"])
        for node in selected_nodes
    }
    failures: list[str] = []
    lockfile = Path(__file__).resolve().parents[1] / "Cargo.lock"
    lockfile_sha256 = canonical_text_sha256(lockfile)
    if lockfile_sha256 != AUDITED_LOCKFILE_SHA256:
        failures.append(
            "lockfile identity is not in the audited pure-Rust envelope: "
            f"{lockfile_sha256}"
        )

    root = package_by_id[metadata["resolve"]["root"]]
    if root["name"] != "georbf":
        failures.append(f"unexpected product package: {root['name']}")
    if root["edition"] != "2024" or root["rust_version"] != "1.85":
        failures.append(
            "product manifest must fix edition 2024 and rust-version 1.85"
        )
    if root["features"]:
        failures.append(f"public Cargo features are forbidden: {root['features']}")

    direct_dependencies = {
        dependency["name"] for dependency in root["dependencies"] if dependency["kind"] is None
    }
    if direct_dependencies != {"faer"}:
        failures.append(f"unexpected production dependencies: {sorted(direct_dependencies)}")

    for package in selected_packages:
        name = package["name"]
        if name in FORBIDDEN_PACKAGES:
            failures.append(f"forbidden native/toolchain package selected: {name}")
        if package.get("links"):
            failures.append(f"package declares native links={package['links']}: {name}")
        license_expression = package.get("license")
        if name != "georbf" and license_expression not in ALLOWED_LICENSE_EXPRESSIONS:
            failures.append(f"unapproved license expression for {name}: {license_expression}")

    for node in selected_nodes:
        package = package_by_id[node["id"]]
        for feature in node["features"]:
            if any(fragment in feature.lower() for fragment in FORBIDDEN_FEATURE_FRAGMENTS):
                failures.append(f"forbidden feature selected: {package['name']}/{feature}")

    if selected_features.get("faer") != ["linalg", "std"]:
        failures.append(f"unexpected faer features: {selected_features.get('faer')}")
    selected_faer = [
        package for package in selected_packages if package["name"] == "faer"
    ]
    if len(selected_faer) != 1 or selected_faer[0]["version"] != "0.24.4":
        failures.append(
            "unexpected faer version: "
            + ",".join(package["version"] for package in selected_faer)
        )

    build_scripts = sorted(
        f"{package['name']} {package['version']}"
        for package in selected_packages
        if any("custom-build" in target["kind"] for target in package["targets"])
    )
    unaudited_build_scripts = set(build_scripts) - AUDITED_BUILD_SCRIPTS
    for build_script in sorted(unaudited_build_scripts):
        failures.append(f"unaudited build script selected: {build_script}")
    license_expressions = sorted(
        {
            package["license"]
            for package in selected_packages
            if package["name"] != "georbf"
        }
    )

    if failures:
        print("audit.result=FAILED")
        for failure in failures:
            print(f"audit.failure={failure}")
        return 1

    print("audit.result=PROVEN")
    print(f"audit.rustc={rustc_version}")
    print(f"audit.target={host}")
    print(f"audit.packages={len(selected_packages)}")
    print(f"audit.lockfile.sha256={lockfile_sha256}")
    print("audit.product.dependencies=faer")
    print(f"audit.faer.version={selected_faer[0]['version']}")
    print("audit.faer.features=linalg,std")
    print("audit.native_links=none")
    print("audit.forbidden_native_packages=none")
    print("audit.forbidden_features=none")
    print(f"audit.build_scripts={';'.join(build_scripts)}")
    print(f"audit.licenses={';'.join(license_expressions)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
