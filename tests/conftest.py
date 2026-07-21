"""pytest configuration: make resources/home/dnanexus importable as modules."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "resources" / "home" / "dnanexus"))
