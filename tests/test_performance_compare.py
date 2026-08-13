"""Tests for the deterministic performance result parser and gate."""

from __future__ import annotations

import importlib.util
import pathlib
import sys
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "compare_performance.py"


def load_comparator():
    spec = importlib.util.spec_from_file_location("compare_performance", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load comparator from {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def output(implementation: str, nanoseconds: int) -> str:
    lines = [
        "georbf-performance-v1 "
        f"implementation={implementation} case=fixed fixed_multi_threads=4 "
        "samples=3 warmups=1 dataset_checksum=1234",
        f"evidence implementation={implementation} "
        "scalars=1.0,2.0,3.0 gradients=1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0",
    ]
    for threads in (1, 4):
        for stage in (
            "preprocess",
            "assembly",
            "solve",
            "scalar_evaluation",
            "gradient_evaluation",
            "end_to_end",
        ):
            for index in range(3):
                lines.append(
                    f"sample implementation={implementation} case=fixed threads={threads} "
                    f"stage={stage} index={index} nanoseconds={nanoseconds + index} "
                    "checksum=abcd"
                )
    return "\n".join(lines)


class PerformanceCompareTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.comparator = load_comparator()

    def test_all_groups_at_or_below_surfe_pass(self) -> None:
        surfe = self.comparator.parse_results(output("surfe", 100))
        georbf = self.comparator.parse_results(output("georbf", 99))
        rows, passed = self.comparator.compare(surfe, georbf)
        self.assertTrue(passed)
        self.assertEqual(len(rows), 14)

    def test_one_slower_group_fails(self) -> None:
        surfe = self.comparator.parse_results(output("surfe", 100))
        georbf = self.comparator.parse_results(output("georbf", 102))
        _, passed = self.comparator.compare(surfe, georbf)
        self.assertFalse(passed)

    def test_rounds_merge_all_samples_and_keep_evidence(self) -> None:
        first = self.comparator.parse_results(output("surfe", 100))
        second = self.comparator.parse_results(output("surfe", 110))
        merged = self.comparator.merge_results([first, second])
        self.assertEqual(merged.header["samples"], "6")
        self.assertEqual(len(merged.samples[(1, "solve")]), 6)
        self.assertEqual(merged.evidence, first.evidence)

    def test_unstable_checksum_is_rejected(self) -> None:
        broken = output("surfe", 100).replace("checksum=abcd", "checksum=bad", 1)
        with self.assertRaisesRegex(ValueError, "unstable checksum"):
            self.comparator.parse_results(broken)


if __name__ == "__main__":
    unittest.main()
