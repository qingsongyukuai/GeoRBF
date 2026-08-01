#!/usr/bin/env python3
"""Generate deterministic high-precision fixtures from declarative cases."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import tempfile
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
CASES_DIR = ROOT / "cases"
FIXTURES_DIR = ROOT / "fixtures"
MANIFEST_PATH = ROOT / "manifest.json"

CASE_SCHEMA_VERSION = "georbf-oracle-case-v1"
FIXTURE_SCHEMA_VERSION = "georbf-oracle-fixture-v1"
MANIFEST_SCHEMA_VERSION = "georbf-oracle-manifest-v1"
GENERATOR_VERSION = "risk-spike-15-v1"
PYTHON_VERSION = "3.12.3"
WORKING_PRECISION_DIGITS = 120
OUTPUT_SIGNIFICANT_DIGITS = 110
ROUNDING_MODE = "ROUND_HALF_EVEN"
OCI_IMAGE = (
    "python:3.12.3-slim-bookworm@"
    "sha256:afc139a0a640942491ec481ad8dda10f2c5b753f5c969393b12480155fe15a63"
)


class GenerationError(RuntimeError):
    pass


def canonical_json_bytes(value: Any, *, pretty: bool = False) -> bytes:
    options: dict[str, Any] = {
        "ensure_ascii": False,
        "sort_keys": True,
    }
    if pretty:
        options["indent"] = 2
    else:
        options["separators"] = (",", ":")
    return (json.dumps(value, **options) + "\n").encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return "sha256:" + hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def as_decimal(value: Any, role: str) -> Decimal:
    if not isinstance(value, str):
        raise GenerationError(f"{role} must be a decimal string")
    try:
        result = Decimal(value)
    except Exception as error:
        raise GenerationError(f"{role} is not a decimal: {value!r}") from error
    if not result.is_finite():
        raise GenerationError(f"{role} must be finite")
    return result


def vector(value: Any, role: str) -> list[Decimal]:
    if not isinstance(value, list) or len(value) != 3:
        raise GenerationError(f"{role} must contain exactly three coordinates")
    return [as_decimal(component, f"{role}[{index}]") for index, component in enumerate(value)]


def matrix(value: Any, role: str) -> list[list[Decimal]]:
    if not isinstance(value, list) or len(value) != 3:
        raise GenerationError(f"{role} must contain exactly three rows")
    return [vector(row, f"{role}[{index}]") for index, row in enumerate(value)]


def dot(left: list[Decimal], right: list[Decimal]) -> Decimal:
    return sum((a * b for a, b in zip(left, right, strict=True)), Decimal(0))


def matvec(coefficients: list[list[Decimal]], operand: list[Decimal]) -> list[Decimal]:
    return [dot(row, operand) for row in coefficients]


def quadratic_form(operand: list[Decimal], coefficients: list[list[Decimal]]) -> Decimal:
    return dot(operand, matvec(coefficients, operand))


def cubic_jet(
    x: list[Decimal], y: list[Decimal], metric: list[list[Decimal]]
) -> dict[str, Any]:
    delta = [a - b for a, b in zip(x, y, strict=True)]
    metric_delta = matvec(metric, delta)
    radius_squared = quadratic_form(delta, metric)
    if radius_squared < 0:
        raise GenerationError("metric produced a negative squared radius")

    if radius_squared.is_zero():
        radius = Decimal(0)
        value = Decimal(0)
        gradient_x = [Decimal(0)] * 3
        mixed_xy = [[Decimal(0)] * 3 for _ in range(3)]
    else:
        radius = radius_squared.sqrt()
        value = radius_squared * radius
        gradient_x = [Decimal(3) * radius * component for component in metric_delta]
        mixed_xy = [
            [
                -Decimal(3) * radius * metric[row][column]
                - Decimal(3) * metric_delta[row] * metric_delta[column] / radius
                for column in range(3)
            ]
            for row in range(3)
        ]

    return {
        "delta": delta,
        "metric_delta": metric_delta,
        "mixed_xy": mixed_xy,
        "radius": radius,
        "radius_squared": radius_squared,
        "value": value,
        "gradient_x": gradient_x,
        "gradient_y": [-component for component in gradient_x],
    }


def evaluate_affine_field(
    field: dict[str, Any],
    support: list[Decimal],
    value_coefficient: Decimal,
    gradient_coefficient: list[Decimal],
) -> Decimal:
    constant = as_decimal(field["constant"], "manufactured_affine_field.constant")
    gradient = vector(field["gradient"], "manufactured_affine_field.gradient")
    field_value = constant + dot(gradient, support)
    return value_coefficient * field_value + dot(gradient_coefficient, gradient)


def generalized_functional(case_input: dict[str, Any]) -> dict[str, Any]:
    metric = matrix(case_input["metric"], "input.metric")
    left = case_input["left"]
    right = case_input["right"]
    left_support = vector(left["support"], "input.left.support")
    right_support = vector(right["support"], "input.right.support")
    left_value_coefficient = as_decimal(
        left["value_coefficient"], "input.left.value_coefficient"
    )
    right_value_coefficient = as_decimal(
        right["value_coefficient"], "input.right.value_coefficient"
    )
    left_gradient_coefficient = vector(
        left["gradient_coefficient"], "input.left.gradient_coefficient"
    )
    right_gradient_coefficient = vector(
        right["gradient_coefficient"], "input.right.gradient_coefficient"
    )
    jet = cubic_jet(left_support, right_support, metric)

    value_value = left_value_coefficient * right_value_coefficient * jet["value"]
    derivative_value = right_value_coefficient * dot(
        left_gradient_coefficient, jet["gradient_x"]
    )
    value_derivative = left_value_coefficient * dot(
        right_gradient_coefficient, jet["gradient_y"]
    )
    derivative_derivative = sum(
        (
            left_gradient_coefficient[row]
            * jet["mixed_xy"][row][column]
            * right_gradient_coefficient[column]
            for row in range(3)
            for column in range(3)
        ),
        Decimal(0),
    )
    pairing = value_value + derivative_value + value_derivative + derivative_derivative
    manufactured_field = case_input["manufactured_affine_field"]

    return {
        "cubic_pairing": pairing,
        "pairing_contributions": {
            "derivative_derivative": derivative_derivative,
            "derivative_value": derivative_value,
            "value_derivative": value_derivative,
            "value_value": value_value,
        },
        "manufactured_affine_observations": {
            "left": evaluate_affine_field(
                manufactured_field,
                left_support,
                left_value_coefficient,
                left_gradient_coefficient,
            ),
            "right": evaluate_affine_field(
                manufactured_field,
                right_support,
                right_value_coefficient,
                right_gradient_coefficient,
            ),
        },
    }


def encode_scalar(value: Decimal) -> dict[str, str]:
    if value.is_zero():
        value = Decimal(0)
    decimal_text = format(value, f".{OUTPUT_SIGNIFICANT_DIGITS - 1}E")
    as_float = 0.0 if value.is_zero() else float(value)
    return {
        "decimal": decimal_text,
        "f64_hex": as_float.hex(),
    }


def encode_result(value: Any) -> Any:
    if isinstance(value, Decimal):
        return encode_scalar(value)
    if isinstance(value, list):
        return [encode_result(item) for item in value]
    if isinstance(value, dict):
        return {key: encode_result(item) for key, item in value.items()}
    raise GenerationError(f"unsupported result type: {type(value).__name__}")


def calculate(case: dict[str, Any]) -> dict[str, Any]:
    case_input = case["input"]
    kind = case["kind"]
    if kind in {"cubic_general_jet", "cubic_origin"}:
        result = cubic_jet(
            vector(case_input["x"], "input.x"),
            vector(case_input["y"], "input.y"),
            matrix(case_input["metric"], "input.metric"),
        )
    elif kind == "cubic_generalized_functional":
        result = generalized_functional(case_input)
    else:
        raise GenerationError(f"unsupported case kind: {kind!r}")
    return encode_result(result)


def load_case(path: Path) -> dict[str, Any]:
    try:
        case = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GenerationError(f"cannot read {path.name}: {error}") from error
    if not isinstance(case, dict):
        raise GenerationError(f"{path.name} must contain an object")
    expected_keys = {"input", "kind", "provenance", "schema_version"}
    if set(case) != expected_keys:
        raise GenerationError(
            f"{path.name} keys must be exactly {sorted(expected_keys)}"
        )
    if case["schema_version"] != CASE_SCHEMA_VERSION:
        raise GenerationError(f"{path.name} has an unsupported schema version")
    return case


def make_fixture(case: dict[str, Any]) -> dict[str, Any]:
    case_id = "case-v1-" + sha256_bytes(canonical_json_bytes(case)).removeprefix("sha256:")
    result = calculate(case)
    payload = {
        "case_id": case_id,
        "input": case["input"],
        "kind": case["kind"],
        "precision": {
            "output_significant_digits": OUTPUT_SIGNIFICANT_DIGITS,
            "rounding": ROUNDING_MODE,
            "working_decimal_digits": WORKING_PRECISION_DIGITS,
        },
        "provenance": case["provenance"],
        "result": result,
        "schema_version": FIXTURE_SCHEMA_VERSION,
    }
    payload["output_sha256"] = sha256_bytes(canonical_json_bytes(result))
    payload["content_sha256"] = sha256_bytes(canonical_json_bytes(payload))
    return payload


def generate(output_root: Path) -> None:
    assert_runtime()
    output_fixtures = output_root / "fixtures"
    output_fixtures.mkdir(parents=True, exist_ok=True)
    case_entries: list[dict[str, Any]] = []

    with localcontext() as context:
        context.prec = WORKING_PRECISION_DIGITS
        context.rounding = ROUND_HALF_EVEN
        for case_path in sorted(CASES_DIR.glob("*.json")):
            case = load_case(case_path)
            fixture = make_fixture(case)
            fixture_name = case_path.name
            fixture_path = output_fixtures / fixture_name
            fixture_bytes = canonical_json_bytes(fixture, pretty=True)
            fixture_path.write_bytes(fixture_bytes)
            case_entries.append(
                {
                    "case_id": fixture["case_id"],
                    "case_path": f"cases/{case_path.name}",
                    "case_sha256": sha256_file(case_path),
                    "fixture_path": f"fixtures/{fixture_name}",
                    "fixture_sha256": sha256_bytes(fixture_bytes),
                    "kind": fixture["kind"],
                    "output_sha256": fixture["output_sha256"],
                    "provenance": fixture["provenance"],
                }
            )

    if not case_entries:
        raise GenerationError("no declarative cases found")

    case_entries.sort(key=lambda entry: entry["case_id"])
    manifest = {
        "cases": case_entries,
        "dependencies": {
            "lockfile": "requirements.lock",
            "lockfile_sha256": sha256_file(ROOT / "requirements.lock"),
            "packages": [],
        },
        "generator": {
            "path": "generate.py",
            "sha256": sha256_file(ROOT / "generate.py"),
            "version": GENERATOR_VERSION,
        },
        "interpreter": {
            "implementation": "CPython",
            "version": PYTHON_VERSION,
        },
        "oci_image": OCI_IMAGE,
        "precision": {
            "output_significant_digits": OUTPUT_SIGNIFICANT_DIGITS,
            "rounding": ROUNDING_MODE,
            "working_decimal_digits": WORKING_PRECISION_DIGITS,
        },
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "verifier": {
            "path": "verify.py",
            "sha256": sha256_file(ROOT / "verify.py"),
        },
    }
    (output_root / "manifest.json").write_bytes(canonical_json_bytes(manifest, pretty=True))


def assert_runtime() -> None:
    actual = platform.python_version()
    if platform.python_implementation() != "CPython" or actual != PYTHON_VERSION:
        raise GenerationError(
            f"generator requires CPython {PYTHON_VERSION}; found "
            f"{platform.python_implementation()} {actual}"
        )


def check_committed_outputs() -> None:
    with tempfile.TemporaryDirectory(prefix="georbf-oracle-") as temporary:
        regenerated_root = Path(temporary)
        generate(regenerated_root)
        expected_paths = [MANIFEST_PATH, *sorted(FIXTURES_DIR.glob("*.json"))]
        actual_paths = [
            regenerated_root / "manifest.json",
            *sorted((regenerated_root / "fixtures").glob("*.json")),
        ]
        expected_names = [path.relative_to(ROOT) for path in expected_paths]
        actual_names = [path.relative_to(regenerated_root) for path in actual_paths]
        if expected_names != actual_names:
            raise GenerationError(
                f"generated output set differs: committed={expected_names}, regenerated={actual_names}"
            )
        changed = [
            str(name)
            for name, expected, actual in zip(
                expected_names, expected_paths, actual_paths, strict=True
            )
            if expected.read_bytes() != actual.read_bytes()
        ]
        if changed:
            raise GenerationError("generated outputs differ: " + ", ".join(changed))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="regenerate in a temporary directory and require byte-for-byte equality",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.check:
            check_committed_outputs()
            print("oracle fixtures regenerate byte-for-byte")
        else:
            generate(ROOT)
            print(f"generated fixtures and {MANIFEST_PATH.name}")
    except GenerationError as error:
        print(f"generation failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
