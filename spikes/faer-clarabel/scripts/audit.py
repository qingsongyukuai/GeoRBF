#!/usr/bin/env python3
"""Audit the active target's locked dependency, feature, build, and license closure."""

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


def run(*command: str) -> str:
    return subprocess.run(command, check=True, text=True, capture_output=True).stdout


def main() -> int:
    host = next(
        line.removeprefix("host: ")
        for line in run("rustc", "--version", "--verbose").splitlines()
        if line.startswith("host: ")
    )
    expected = os.environ.get("EXPECTED_TARGET")
    if expected and expected != host:
        raise RuntimeError(f"runner target mismatch: expected {expected}, found {host}")

    cargo_command = ["cargo"]
    sparse_mirror = os.environ.get("GEORBF_CARGO_SPARSE_REGISTRY")
    if sparse_mirror:
        cargo_command.extend(
            [
                "--config",
                'source.crates-io.replace-with="georbf-probe-mirror"',
                "--config",
                f'source.georbf-probe-mirror.registry="sparse+{sparse_mirror.rstrip("/")}/"',
            ]
        )

    metadata = json.loads(
        run(
            *cargo_command,
            "metadata",
            "--locked",
            "--format-version",
            "1",
            "--filter-platform",
            host,
        )
    )
    packages = metadata["packages"]
    package_by_id = {package["id"]: package for package in packages}
    failures: list[str] = []

    for package in packages:
        name = package["name"]
        if name in FORBIDDEN_PACKAGES:
            failures.append(f"forbidden native/toolchain package selected: {name}")
        if package.get("links"):
            failures.append(f"package declares native links={package['links']}: {name}")
        license_expression = package.get("license")
        if name != "georbf-backend-probe" and license_expression not in ALLOWED_LICENSE_EXPRESSIONS:
            failures.append(f"unapproved license expression for {name}: {license_expression}")

    for node in metadata["resolve"]["nodes"]:
        package = package_by_id[node["id"]]
        for feature in node["features"]:
            if any(fragment in feature.lower() for fragment in FORBIDDEN_FEATURE_FRAGMENTS):
                failures.append(f"forbidden feature selected: {package['name']}/{feature}")

    selected = {
        package_by_id[node["id"]]["name"]: node["features"]
        for node in metadata["resolve"]["nodes"]
    }
    if selected.get("faer") != ["linalg", "std"]:
        failures.append(f"unexpected faer features: {selected.get('faer')}")
    if selected.get("clarabel") != ["serde"]:
        failures.append(f"unexpected Clarabel features: {selected.get('clarabel')}")

    build_scripts = sorted(
        f"{package['name']} {package['version']}"
        for package in packages
        if any("custom-build" in target["kind"] for target in package["targets"])
    )
    license_expressions = sorted(
        {
            package["license"]
            for package in packages
            if package["name"] != "georbf-backend-probe"
        }
    )

    if failures:
        print("audit.result=FAILED")
        for failure in failures:
            print(f"audit.failure={failure}")
        return 1

    print("audit.result=PROVEN")
    print(f"audit.target={host}")
    print(f"audit.packages={len(packages)}")
    lockfile = Path(__file__).resolve().parents[1] / "Cargo.lock"
    print(f"audit.lockfile.sha256={hashlib.sha256(lockfile.read_bytes()).hexdigest()}")
    print("audit.native_links=none")
    print("audit.forbidden_native_packages=none")
    print("audit.forbidden_features=none")
    print(f"audit.build_scripts={';'.join(build_scripts)}")
    print(f"audit.licenses={';'.join(license_expressions)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
