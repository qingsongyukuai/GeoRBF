#!/usr/bin/env python3
"""Generate issue #41's independent precision-rescue oracle corpus."""

from __future__ import annotations

import hashlib
import json
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CASE_PATH = ROOT / "cases" / "precision-rescue.json"
FIXTURE_PATH = ROOT / "fixtures" / "precision-rescue.json"
MANIFEST_PATH = ROOT / "source-manifest.json"
WORKING_DIGITS = 160
OUTPUT_DIGITS = 120


def sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def decimal(value: str) -> Decimal:
    return Decimal(value)


def dd_value(components: list[str]) -> Decimal:
    return sum((decimal(component) for component in components), Decimal(0))


def decimal_text(value: Decimal) -> str:
    if value.is_zero():
        return "0"
    with localcontext() as context:
        context.prec = OUTPUT_DIGITS
        context.rounding = ROUND_HALF_EVEN
        return format(+value, ".119e")


def scalar(value: Decimal) -> dict[str, str]:
    rounded = float(value)
    return {"decimal": decimal_text(value), "f64_hex": rounded.hex()}


def double_double(value: Decimal) -> dict[str, dict[str, str]]:
    high = float(value)
    residual = value - Decimal.from_float(high)
    low = float(residual)
    return {
        "high": scalar(Decimal.from_float(high)),
        "low": scalar(Decimal.from_float(low)),
        "value": {"decimal": decimal_text(value)},
    }


def arithmetic(case: dict[str, object]) -> dict[str, object]:
    left = dd_value(case["left"])
    operation = case["operation"]
    if operation == "sqrt":
        result = left.sqrt()
    else:
        right = dd_value(case["right"])
        result = {
            "add": left + right,
            "subtract": left - right,
            "multiply": left * right,
            "divide": left / right,
        }[operation]
    return {"operation": operation, "result": double_double(result)}


def dot(left: list[Decimal], right: list[Decimal]) -> Decimal:
    return sum((a * b for a, b in zip(left, right)), Decimal(0))


def matrix_vector(matrix: list[list[Decimal]], vector: list[Decimal]) -> list[Decimal]:
    return [dot(row, vector) for row in matrix]


def cubic_jet(
    left_support: list[Decimal],
    right_support: list[Decimal],
    metric: list[list[Decimal]],
) -> tuple[Decimal, list[Decimal], list[Decimal], list[list[Decimal]]]:
    delta = [left - right for left, right in zip(left_support, right_support)]
    if all(component.is_zero() for component in delta):
        return Decimal(0), [Decimal(0)] * 3, [Decimal(0)] * 3, [[Decimal(0)] * 3 for _ in range(3)]
    metric_delta = matrix_vector(metric, delta)
    radius = dot(delta, metric_delta).sqrt()
    gradient_left = [Decimal(3) * radius * component for component in metric_delta]
    gradient_right = [-component for component in gradient_left]
    mixed = [
        [
            -Decimal(3)
            * (radius * metric[row][column] + metric_delta[row] * metric_delta[column] / radius)
            for column in range(3)
        ]
        for row in range(3)
    ]
    return radius**3, gradient_left, gradient_right, mixed


def term(raw: dict[str, object]) -> tuple[list[Decimal], Decimal, list[Decimal]]:
    return (
        [decimal(value) for value in raw["support"]],
        decimal(raw["value"]),
        [decimal(value) for value in raw["gradient"]],
    )


def cubic_pairing(case: dict[str, object]) -> Decimal:
    metric = [[decimal(value) for value in row] for row in case["metric"]]
    result = Decimal(0)
    for raw_left in case["left"]:
        left_support, left_value, left_gradient = term(raw_left)
        for raw_right in case["right"]:
            right_support, right_value, right_gradient = term(raw_right)
            value, gradient_left_jet, gradient_right_jet, mixed = cubic_jet(
                left_support, right_support, metric
            )
            result += left_value * right_value * value
            result += left_value * dot(right_gradient, gradient_right_jet)
            result += right_value * dot(left_gradient, gradient_left_jet)
            result += dot(left_gradient, matrix_vector(mixed, right_gradient))
    return result


def schur(case: dict[str, object]) -> dict[str, object]:
    value = dd_value(case["diagonal"])
    for left, right in zip(case["left_factors"], case["right_factors"]):
        value -= dd_value(left) * dd_value(right)
    return {"classification": case["classification"], "result": double_double(value)}


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    with localcontext() as context:
        context.prec = WORKING_DIGITS
        context.rounding = ROUND_HALF_EVEN
        declarations = json.loads(CASE_PATH.read_text())
        fixture = {
            "arithmetic": [arithmetic(case) for case in declarations["arithmetic"]],
            "cubic_pairing": double_double(cubic_pairing(declarations["cubic_pairing"])),
            "provenance": {
                "dependency_closure": [],
                "generator": "independent Python decimal arithmetic",
                "rounding": "ROUND_HALF_EVEN",
                "source": "GeoRBF issue #41",
                "working_decimal_digits": WORKING_DIGITS,
            },
            "schema_version": "georbf-precision-rescue-fixture-v1",
            "schur": [schur(case) for case in declarations["schur"]],
        }
        FIXTURE_PATH.parent.mkdir(parents=True, exist_ok=True)
        FIXTURE_PATH.write_bytes(canonical_json(fixture))

    manifest = {
        "case": {"path": str(CASE_PATH.relative_to(ROOT)), "sha256": sha256(CASE_PATH)},
        "dependencies": {"packages": []},
        "fixture": {"path": str(FIXTURE_PATH.relative_to(ROOT)), "sha256": sha256(FIXTURE_PATH)},
        "generator": {"path": "generate.py", "sha256": sha256(Path(__file__))},
        "interpreter": {"implementation": "CPython", "version": "3.12.3"},
        "precision": {
            "output_significant_digits": OUTPUT_DIGITS,
            "rounding": "ROUND_HALF_EVEN",
            "working_decimal_digits": WORKING_DIGITS,
        },
        "provenance": {"source": "GeoRBF issue #41"},
        "schema_version": "georbf-precision-rescue-manifest-v1",
    }
    MANIFEST_PATH.write_bytes(canonical_json(manifest))


if __name__ == "__main__":
    main()
