from __future__ import annotations

import json
import hashlib
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return "sha256:" + hashlib.sha256(value).hexdigest()


class OracleReproducibilityTests(unittest.TestCase):
    def run_verifier(self, root: Path = ROOT) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, "verify.py"],
            cwd=root,
            check=False,
            capture_output=True,
            text=True,
        )

    def test_committed_outputs_regenerate_without_a_byte_diff(self) -> None:
        completed = subprocess.run(
            [sys.executable, "generate.py", "--check"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_committed_manifest_hashes_and_f64_encodings_validate(self) -> None:
        completed = self.run_verifier()

        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_worked_examples_have_the_expected_known_f64_values(self) -> None:
        general = json.loads(
            (ROOT / "fixtures/cubic-general-jet.json").read_text(encoding="utf-8")
        )
        functional = json.loads(
            (ROOT / "fixtures/cubic-generalized-functional.json").read_text(
                encoding="utf-8"
            )
        )

        self.assertEqual(
            {
                "cubic_pairing": functional["result"]["cubic_pairing"]["f64_hex"],
                "left_affine_truth": functional["result"][
                    "manufactured_affine_observations"
                ]["left"]["f64_hex"],
                "radius_squared": general["result"]["radius_squared"]["f64_hex"],
                "value": general["result"]["value"]["f64_hex"],
            },
            {
                "cubic_pairing": "-0x1.e8b2ef61dc324p+4",
                "left_affine_truth": "0x1.4000000000000p-1",
                "radius_squared": "0x1.5600000000000p+4",
                "value": "0x1.8b4b0530981ddp+6",
            },
        )

    def test_verifier_stably_rejects_a_tampered_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            copied_root = Path(temporary) / "oracle-fixtures"
            shutil.copytree(ROOT, copied_root)
            fixture_path = copied_root / "fixtures/cubic-general-jet.json"
            fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
            fixture["result"]["value"]["decimal"] = "1.0E+0"
            fixture_path.write_text(
                json.dumps(fixture, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

            completed = self.run_verifier(copied_root)

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("fixture hash mismatch", completed.stderr)

    def test_verifier_rejects_a_short_decimal_even_with_rehashed_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            copied_root = Path(temporary) / "oracle-fixtures"
            shutil.copytree(ROOT, copied_root)
            fixture_path = copied_root / "fixtures/cubic-general-jet.json"
            fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
            fixture["result"]["value"]["decimal"] = "9.8823261985096404E+1"
            fixture["output_sha256"] = sha256_bytes(
                canonical_json_bytes(fixture["result"])
            )
            content = dict(fixture)
            content.pop("content_sha256")
            fixture["content_sha256"] = sha256_bytes(canonical_json_bytes(content))
            fixture_bytes = (
                json.dumps(fixture, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
            ).encode("utf-8")
            fixture_path.write_bytes(fixture_bytes)

            manifest_path = copied_root / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            entry = next(
                item
                for item in manifest["cases"]
                if item["fixture_path"] == "fixtures/cubic-general-jet.json"
            )
            entry["output_sha256"] = fixture["output_sha256"]
            entry["fixture_sha256"] = sha256_bytes(fixture_bytes)
            manifest_path.write_text(
                json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

            completed = self.run_verifier(copied_root)

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("110 significant digits", completed.stderr)

    def test_verifier_rejects_a_changed_oci_pin(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            copied_root = Path(temporary) / "oracle-fixtures"
            shutil.copytree(ROOT, copied_root)
            manifest_path = copied_root / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["oci_image"] = "python:latest"
            manifest_path.write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

            completed = self.run_verifier(copied_root)

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("OCI image pin", completed.stderr)


if __name__ == "__main__":
    unittest.main()
