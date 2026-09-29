import os
import shutil
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

EXAMPLE = os.path.join(ROOT, "examples", "minimal", "bundle")
DATA = os.path.join(os.path.dirname(__file__), "data")
SCRIPTS = os.path.join(ROOT, "scripts")

needs_htslib = pytest.mark.skipif(not (shutil.which("bgzip") and shutil.which("tabix")),
                                  reason="bgzip/tabix (htslib) not on PATH: run `module load htslib` or install htslib")


@pytest.fixture(scope="session")
def example_bundle():
    return EXAMPLE


@pytest.fixture
def tmp_bundle(tmp_path):
    """A writable copy of the example bundle."""
    dst = tmp_path / "bundle"
    shutil.copytree(EXAMPLE, dst)
    return str(dst)
