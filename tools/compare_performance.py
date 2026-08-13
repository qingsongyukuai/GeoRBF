#!/usr/bin/env python3
"""Run and compare the frozen Surfe and GeoRBF performance harnesses."""

from __future__ import annotations

import argparse
import math
import shlex
import statistics
import subprocess
import sys
from dataclasses import dataclass


STAGES = (
    "preprocess",
    "assembly",
    "solve",
    "scalar_evaluation",
    "gradient_evaluation",
    "end_to_end",
)
THREADED_STAGES = {"scalar_evaluation", "gradient_evaluation", "end_to_end"}


@dataclass(frozen=True)
class Results:
    header: dict[str, str]
    samples: dict[tuple[int, str], tuple[int, ...]]
    checksums: dict[tuple[int, str], str]
    evidence: dict[str, tuple[float, ...]]


def fields(line: str) -> dict[str, str]:
    return dict(token.split("=", 1) for token in line.split()[1:] if "=" in token)


def parse_results(output: str) -> Results:
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    header_line = next(
        (line for line in lines if line.startswith("georbf-performance-v1 ")), None
    )
    if header_line is None:
        raise ValueError("benchmark header is missing")
    header = fields(header_line)
    evidence_line = next((line for line in lines if line.startswith("evidence ")), None)
    if evidence_line is None:
        raise ValueError("benchmark evidence is missing")
    evidence_fields = fields(evidence_line)
    evidence = {
        name: tuple(float(value) for value in evidence_fields[name].split(","))
        for name in ("scalars", "gradients")
    }
    if len(evidence["scalars"]) != 3 or len(evidence["gradients"]) != 9:
        raise ValueError("benchmark evidence has the wrong shape")
    grouped: dict[tuple[int, str], list[int]] = {}
    checksums: dict[tuple[int, str], str] = {}
    for line in lines:
        if not line.startswith("sample "):
            continue
        row = fields(line)
        key = (int(row["threads"]), row["stage"])
        grouped.setdefault(key, []).append(int(row["nanoseconds"]))
        previous = checksums.setdefault(key, row["checksum"])
        if previous != row["checksum"]:
            raise ValueError(f"unstable checksum for threads={key[0]} stage={key[1]}")
    expected_samples = int(header["samples"])
    fixed_multi_threads = int(header["fixed_multi_threads"])
    expected_keys = {
        (threads, stage) for threads in (1, fixed_multi_threads) for stage in STAGES
    }
    if set(grouped) != expected_keys:
        raise ValueError("benchmark did not emit all 12 thread/stage groups")
    for key, values in grouped.items():
        if len(values) != expected_samples:
            raise ValueError(f"wrong sample count for threads={key[0]} stage={key[1]}")
    return Results(
        header,
        {key: tuple(values) for key, values in grouped.items()},
        checksums,
        evidence,
    )


def run(command: str) -> Results:
    completed = subprocess.run(
        shlex.split(command),
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.stderr:
        sys.stderr.write(completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(f"benchmark failed ({completed.returncode}): {command}")
    return parse_results(completed.stdout)


def merge_results(results: list[Results]) -> Results:
    if not results:
        raise ValueError("at least one benchmark round is required")
    first = results[0]
    for result in results[1:]:
        for key in ("implementation", "case", "fixed_multi_threads", "warmups", "dataset_checksum"):
            if result.header.get(key) != first.header.get(key):
                raise ValueError(f"benchmark round mismatch for {key}")
        if result.checksums != first.checksums or result.evidence != first.evidence:
            raise ValueError("benchmark round evidence is unstable")
    header = dict(first.header)
    header["samples"] = str(sum(int(result.header["samples"]) for result in results))
    samples = {
        key: tuple(value for result in results for value in result.samples[key])
        for key in first.samples
    }
    return Results(header, samples, first.checksums, first.evidence)


def compare(surfe: Results, georbf: Results) -> tuple[list[str], bool]:
    for key in ("case", "fixed_multi_threads", "samples", "warmups", "dataset_checksum"):
        if surfe.header.get(key) != georbf.header.get(key):
            raise ValueError(f"benchmark header mismatch for {key}")
    if surfe.checksums[(1, "assembly")] != georbf.checksums[(1, "assembly")]:
        raise ValueError("assembled matrix checksum mismatch")
    for name, absolute, relative in (
        ("scalars", 1.0e-9, 1.0e-8),
        ("gradients", 1.0e-8, 1.0e-7),
    ):
        for index, (expected, actual) in enumerate(
            zip(surfe.evidence[name], georbf.evidence[name])
        ):
            if not math.isclose(actual, expected, abs_tol=absolute, rel_tol=relative):
                raise ValueError(f"{name} parity mismatch at index {index}")

    rows = [
        "| Threads | Stage | Surfe median (ns) | GeoRBF median (ns) | Ratio | Gate |",
        "|---:|---|---:|---:|---:|:---:|",
    ]
    passed = True
    fixed_multi_threads = int(surfe.header["fixed_multi_threads"])
    for threads in (1, fixed_multi_threads):
        for stage in STAGES:
            measurement_threads = threads if stage in THREADED_STAGES else 1
            key = (measurement_threads, stage)
            surfe_median = int(statistics.median(surfe.samples[key]))
            georbf_median = int(statistics.median(georbf.samples[key]))
            gate = georbf_median <= surfe_median
            passed &= gate
            rows.append(
                f"| {threads} | {stage} | {surfe_median} | {georbf_median} | "
                f"{georbf_median / surfe_median:.3f} | {'PASS' if gate else 'FAIL'} |"
            )
    return rows, passed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--surfe-command", required=True)
    parser.add_argument("--georbf-command", required=True)
    parser.add_argument("--rounds", type=int, default=3)
    arguments = parser.parse_args()
    if arguments.rounds < 1:
        parser.error("--rounds must be positive")
    surfe_rounds = []
    georbf_rounds = []
    for round_index in range(arguments.rounds):
        if round_index % 2 == 0:
            surfe_rounds.append(run(arguments.surfe_command))
            georbf_rounds.append(run(arguments.georbf_command))
        else:
            georbf_rounds.append(run(arguments.georbf_command))
            surfe_rounds.append(run(arguments.surfe_command))
    rows, passed = compare(
        merge_results(surfe_rounds),
        merge_results(georbf_rounds),
    )
    print("\n".join(rows))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
