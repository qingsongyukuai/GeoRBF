import hashlib
import tempfile
import unittest
from pathlib import Path

from audit import byte_sha256


class ByteSha256Tests(unittest.TestCase):
    def test_line_endings_are_part_of_the_identity(self) -> None:
        lf_content = b"version = 4\n\n[[package]]\nname = \"example\"\n"
        expected = hashlib.sha256(lf_content).hexdigest()

        with tempfile.TemporaryDirectory() as directory:
            lockfile = Path(directory) / "Cargo.lock"
            lockfile.write_bytes(lf_content)
            self.assertEqual(byte_sha256(lockfile), expected)

            lockfile.write_bytes(lf_content.replace(b"\n", b"\r\n"))
            self.assertNotEqual(byte_sha256(lockfile), expected)


if __name__ == "__main__":
    unittest.main()
