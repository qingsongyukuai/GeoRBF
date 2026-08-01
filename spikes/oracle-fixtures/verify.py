#!/usr/bin/env python3
"""Verify fixture identities, hashes, manifest coverage, and exact f64 encodings."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any


DEFAULT_ROOT = Path(__file__).resolve().parent
CASE_SCHEMA_VERSION = "georbf-oracle-case-v1"
FIXTURE_SCHEMA_VERSION = "georbf-oracle-fixture-v1"
MANIFEST_SCHEMA_VERSION = "georbf-oracle-manifest-v1"
EXPECTED_INTERPRETER = {"implementation": "CPython", "version": "3.12.3"}
EXPECTED_OCI_IMAGE = (
    "python:3.12.3-slim-bookworm@"
    "sha256:afc139a0a640942491ec481ad8dda10f2c5b753f5c969393b12480155fe15a63"
)
EXPECTED_PRECISION = {
    "output_significant_digits": 110,
    "rounding": "ROUND_HALF_EVEN",
    "working_decimal_digits": 120,
}
EXPECTED_GENERATOR_VERSION = "risk-spike-15-v1"
DECIMAL_ENCODING = re.compile(r"-?\d\.\d{109}E[+-]\d+")


class VerificationError(RuntimeError):
    pass


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return "sha256:" + hashlib.sha256(value).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise VerificationError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise VerificationError(f"{path} must contain an object")
    return value


def verify_encoded_scalars(value: Any, role: str) -> None:
    if isinstance(value, list):
        for index, item in enumerate(value):
            verify_encoded_scalars(item, f"{role}[{index}]")
        return
    if not isinstance(value, dict):
        raise VerificationError(f"{role} is not an encoded scalar or collection")
    if set(value) == {"decimal", "f64_hex"}:
        if not isinstance(value["decimal"], str) or DECIMAL_ENCODING.fullmatch(
            value["decimal"]
        ) is None:
            raise VerificationError(f"{role} must contain exactly 110 significant digits")
        try:
            decimal_value = Decimal(value["decimal"])
            hex_value = float.fromhex(value["f64_hex"])
        except (ArithmeticError, TypeError, ValueError) as error:
            raise VerificationError(f"{role} has an invalid numeric encoding") from error
        rounded = 0.0 if decimal_value.is_zero() else float(decimal_value)
        if not decimal_value.is_finite() or not math.isfinite(hex_value):
            raise VerificationError(f"{role} must be finite")
        if rounded.hex() != hex_value.hex() or hex_value.hex() != value["f64_hex"]:
            raise VerificationError(f"{role} decimal and f64 hex do not round-trip")
        return
    for key, item in value.items():
        verify_encoded_scalars(item, f"{role}.{key}")


def verify_fixture(root: Path, entry: dict[str, Any]) -> None:
    case_path = root / entry["case_path"]
    fixture_path = root / entry["fixture_path"]
    if sha256_bytes(case_path.read_bytes()) != entry["case_sha256"]:
        raise VerificationError(f"case hash mismatch: {entry['case_path']}")
    if sha256_bytes(fixture_path.read_bytes()) != entry["fixture_sha256"]:
        raise VerificationError(f"fixture hash mismatch: {entry['fixture_path']}")

    case = load_json(case_path)
    fixture = load_json(fixture_path)
    if case.get("schema_version") != CASE_SCHEMA_VERSION:
        raise VerificationError(f"unsupported case schema: {entry['case_path']}")
    if fixture.get("schema_version") != FIXTURE_SCHEMA_VERSION:
        raise VerificationError(f"unsupported fixture schema: {entry['fixture_path']}")
    if fixture.get("precision") != EXPECTED_PRECISION:
        raise VerificationError(f"fixture precision mismatch: {entry['fixture_path']}")
    expected_case_id = "case-v1-" + sha256_bytes(
        canonical_json_bytes(case)
    ).removeprefix("sha256:")
    if fixture.get("case_id") != expected_case_id or entry["case_id"] != expected_case_id:
        raise VerificationError(f"CaseId mismatch: {entry['fixture_path']}")
    if fixture.get("input") != case.get("input"):
        raise VerificationError(f"canonical input mismatch: {entry['fixture_path']}")
    if fixture.get("kind") != case.get("kind") or entry["kind"] != case.get("kind"):
        raise VerificationError(f"case kind mismatch: {entry['fixture_path']}")
    if fixture.get("provenance") != case.get("provenance"):
        raise VerificationError(f"provenance mismatch: {entry['fixture_path']}")

    result = fixture.get("result")
    output_hash = sha256_bytes(canonical_json_bytes(result))
    if fixture.get("output_sha256") != output_hash or entry["output_sha256"] != output_hash:
        raise VerificationError(f"output hash mismatch: {entry['fixture_path']}")
    content = dict(fixture)
    content_hash = content.pop("content_sha256", None)
    if content_hash != sha256_bytes(canonical_json_bytes(content)):
        raise VerificationError(f"content hash mismatch: {entry['fixture_path']}")
    verify_encoded_scalars(result, f"fixture[{entry['case_id']}].result")


def verify(root: Path) -> None:
    manifest_path = root / "manifest.json"
    manifest = load_json(manifest_path)
    if manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        raise VerificationError("unsupported manifest schema")
    if manifest.get("interpreter") != EXPECTED_INTERPRETER:
        raise VerificationError("interpreter pin mismatch")
    if manifest.get("oci_image") != EXPECTED_OCI_IMAGE:
        raise VerificationError("OCI image pin mismatch")
    if manifest.get("precision") != EXPECTED_PRECISION:
        raise VerificationError("precision policy mismatch")
    if manifest.get("dependencies", {}).get("packages") != []:
        raise VerificationError("dependency package set mismatch")
    if manifest.get("generator", {}).get("version") != EXPECTED_GENERATOR_VERSION:
        raise VerificationError("generator version mismatch")
    for section in ("generator", "verifier"):
        source = root / manifest[section]["path"]
        if sha256_bytes(source.read_bytes()) != manifest[section]["sha256"]:
            raise VerificationError(f"{section} source hash mismatch")
    lock = root / manifest["dependencies"]["lockfile"]
    if sha256_bytes(lock.read_bytes()) != manifest["dependencies"]["lockfile_sha256"]:
        raise VerificationError("dependency lockfile hash mismatch")
    entries = manifest.get("cases")
    if not isinstance(entries, list) or not entries:
        raise VerificationError("manifest has no cases")
    case_ids = [entry["case_id"] for entry in entries]
    if case_ids != sorted(case_ids) or len(case_ids) != len(set(case_ids)):
        raise VerificationError("manifest CaseIds are not sorted and unique")
    declared = sorted(path.relative_to(root).as_posix() for path in (root / "cases").glob("*.json"))
    fixed = sorted(path.relative_to(root).as_posix() for path in (root / "fixtures").glob("*.json"))
    if declared != sorted(entry["case_path"] for entry in entries):
        raise VerificationError("manifest does not cover the declarative cases exactly")
    if fixed != sorted(entry["fixture_path"] for entry in entries):
        raise VerificationError("manifest does not cover the fixtures exactly")
    for entry in entries:
        verify_fixture(root, entry)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="oracle spike root containing manifest.json",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        verify(args.root.resolve())
        print("oracle manifest, hashes, and f64 encodings are valid")
    except (KeyError, OSError, VerificationError) as error:
        print(f"verification failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
