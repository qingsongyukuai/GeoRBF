"""Validate one requirement-level trace for every v0.2.0 contract."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


TRACEABILITY_PATH = Path("validation/v0.2.0/traceability.json")
REFERENCE_SET_KINDS = ("evidence_sets", "documentation_sets")


def requirement_ids(prefix: str, numbers: list[int]) -> set[str]:
    return {f"{prefix}-{number:03}" for number in numbers}


EXPECTED_REQUIREMENTS = frozenset(
    requirement_ids("PAPI", [*range(1, 16), 17, 18, 19])
    | requirement_ids("DOM", [*range(1, 8), *range(9, 23)])
    | requirement_ids("IR", [*range(1, 15)])
    | requirement_ids("KER", [1, *range(4, 10)])
    | requirement_ids("NUM", [*range(1, 16)])
    | requirement_ids("DIA", [*range(1, 10)])
    | requirement_ids("VAL", [*range(1, 16)])
)


def safe_repository_path(root: Path, value: object) -> tuple[Path | None, str | None]:
    if not isinstance(value, str) or not value:
        return None, "reference path must be a non-empty string"
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        return None, f"reference path escapes the repository: {value}"
    return root / relative, None


def validate_reference(root: Path, role: str, reference: object) -> list[str]:
    if not isinstance(reference, dict):
        return [f"{role} reference must be an object"]
    path, path_failure = safe_repository_path(root, reference.get("path"))
    if path_failure is not None:
        return [f"{role}: {path_failure}"]
    assert path is not None
    if not path.is_file():
        return [f"{role}: missing referenced path {path.relative_to(root).as_posix()}"]
    marker = reference.get("contains")
    if not isinstance(marker, str) or not marker:
        return [f"{role}: marker must be a non-empty string"]
    try:
        content = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        return [f"{role}: cannot read {path.relative_to(root).as_posix()}: {error}"]
    if marker not in content:
        return [
            f"{role}: missing marker {marker!r} in {path.relative_to(root).as_posix()}"
        ]
    return []


def validate_reference_sets(
    root: Path, kind: str, value: object
) -> tuple[set[str], list[str]]:
    failures: list[str] = []
    if not isinstance(value, dict) or not value:
        return set(), [f"traceability {kind} must be a non-empty object"]
    names: set[str] = set()
    for name, references in value.items():
        role = f"{kind}.{name}"
        if not isinstance(name, str) or not name:
            failures.append(f"{kind} contains an empty or non-string set name")
            continue
        names.add(name)
        if not isinstance(references, list) or not references:
            failures.append(f"{role} must be a non-empty list")
            continue
        for reference in references:
            failures.extend(validate_reference(root, role, reference))
    return names, failures


def validate_traceability(root: Path, document: object) -> list[str]:
    failures: list[str] = []
    if not isinstance(document, dict):
        return ["traceability root must be an object"]
    if document.get("schema_version") != "georbf-traceability-v1":
        failures.append("traceability schema_version must be georbf-traceability-v1")
    if document.get("release") != "0.2.0":
        failures.append("traceability release must be 0.2.0")

    reference_names: dict[str, set[str]] = {}
    for kind in REFERENCE_SET_KINDS:
        names, reference_failures = validate_reference_sets(root, kind, document.get(kind))
        reference_names[kind] = names
        failures.extend(reference_failures)

    traces = document.get("requirements")
    if not isinstance(traces, list) or not traces:
        return failures + ["traceability requirements must be a non-empty list"]
    declared_ids: list[str] = []
    behaviors: list[str] = []
    api_ir_paths: list[str] = []
    used_sets = {kind: set() for kind in REFERENCE_SET_KINDS}
    for index, trace in enumerate(traces):
        role = f"requirements[{index}]"
        if not isinstance(trace, dict):
            failures.append(f"{role} must be an object")
            continue
        requirement = trace.get("id")
        if not isinstance(requirement, str):
            failures.append(f"{role}.id must be a string")
        else:
            declared_ids.append(requirement)
        behavior = trace.get("behavior")
        if not isinstance(behavior, str) or not behavior.strip():
            failures.append(f"{role}.behavior must be a non-empty string")
        else:
            behaviors.append(behavior)
        api_ir = trace.get("api_ir")
        if not isinstance(api_ir, str) or not api_ir.strip():
            failures.append(f"{role}.api_ir must be a non-empty public API or IR path")
        else:
            api_ir_paths.append(api_ir)
        for field, kind in (
            ("evidence", "evidence_sets"),
            ("documentation", "documentation_sets"),
        ):
            set_name = trace.get(field)
            if not isinstance(set_name, str) or not set_name:
                failures.append(f"{role}.{field} must name a reference set")
            elif set_name not in reference_names[kind]:
                failures.append(f"{role}.{field} names unknown set {set_name!r}")
            else:
                used_sets[kind].add(set_name)

    counts = Counter(declared_ids)
    duplicates = sorted(requirement for requirement, count in counts.items() if count > 1)
    if duplicates:
        failures.append(f"duplicate requirement traces: {', '.join(duplicates)}")
    declared = set(declared_ids)
    missing = sorted(EXPECTED_REQUIREMENTS - declared)
    unexpected = sorted(declared - EXPECTED_REQUIREMENTS)
    if missing:
        failures.append(f"missing requirements: {', '.join(missing)}")
    if unexpected:
        failures.append(f"out-of-scope requirements: {', '.join(unexpected)}")
    for role, values in (("behaviors", behaviors), ("API/IR paths", api_ir_paths)):
        if len(values) != len(set(values)):
            failures.append(f"requirement-level {role} must be unique")
    for kind in REFERENCE_SET_KINDS:
        unused = sorted(reference_names[kind] - used_sets[kind])
        if unused:
            failures.append(f"unused {kind}: {', '.join(unused)}")
    return failures


def load_traceability(root: Path) -> tuple[object | None, list[str]]:
    path = root / TRACEABILITY_PATH
    try:
        return json.loads(path.read_text(encoding="utf-8")), []
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        return None, [f"cannot read {TRACEABILITY_PATH.as_posix()}: {error}"]
