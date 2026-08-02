import hashlib
import tempfile
import unittest
from pathlib import Path

from audit import canonical_text_sha256


class CanonicalTextSha256Tests(unittest.TestCase):
    def test_line_endings_do_not_change_the_identity(self) -> None:
        lf_content = b"version = 4\n\n[[package]]\nname = \"example\"\n"
        expected = hashlib.sha256(lf_content).hexdigest()

        with tempfile.TemporaryDirectory() as directory:
            lockfile = Path(directory) / "Cargo.lock"
            for content in (lf_content, lf_content.replace(b"\n", b"\r\n")):
                lockfile.write_bytes(content)
                self.assertEqual(canonical_text_sha256(lockfile), expected)


if __name__ == "__main__":
    unittest.main()
