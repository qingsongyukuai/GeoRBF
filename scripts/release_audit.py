#!/usr/bin/env python3
"""Aggregate the fail-closed v0.1.0 Equality Spine release audits."""

from __future__ import annotations

import sys
from pathlib import Path

from release_repository_checks import audit_repository_artifacts
from release_traceability import (
    EXPECTED_REQUIREMENTS,
    load_traceability,
    validate_traceability,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def audit_repository(root: Path) -> list[str]:
    traceability, failures = load_traceability(root)
    if traceability is not None:
        failures.extend(validate_traceability(root, traceability))
    failures.extend(audit_repository_artifacts(root))
    return failures


def main() -> int:
    failures = audit_repository(REPOSITORY_ROOT)
    if failures:
        print("release-audit.result=FAILED")
        for failure in failures:
            print(f"release-audit.failure={failure}")
        return 1
    print("release-audit.result=PROVEN")
    print("release-audit.release=0.1.0")
    print(f"release-audit.requirements={len(EXPECTED_REQUIREMENTS)}")
    print("release-audit.oracle-adoption=byte-identical")
    print("release-audit.placeholders=none")
    return 0


if __name__ == "__main__":
    sys.exit(main())
