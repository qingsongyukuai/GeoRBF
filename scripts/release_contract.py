"""Single release identity consumed by fail-closed repository audits."""

from pathlib import Path


RELEASE_VERSION = "0.2.0"
RELEASE_TAG = f"v{RELEASE_VERSION}"
TRACEABILITY_PATH = Path("validation") / RELEASE_TAG / "traceability.json"
WORKFLOW_PATH = Path(".github/workflows/product-v0.2.yml")
