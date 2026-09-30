import re
import sys
import sysconfig
import tomllib
from pathlib import Path

import pybgzf
import pytest

ROOT = Path(__file__).parent.parent


def test_readme_links_work_on_pypi() -> None:
    targets = re.findall(r"\]\(([^)]+)\)", (ROOT / "README.md").read_text())
    assert targets
    assert all(target.startswith("https://") for target in targets), targets


def test_the_version_comes_from_cargo() -> None:
    cargo = tomllib.loads((ROOT / "Cargo.toml").read_text())
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert "version" not in project
    assert project["dynamic"] == ["version"]
    assert pybgzf.__version__ == cargo["package"]["version"]


@pytest.mark.skipif(
    not sysconfig.get_config_var("Py_GIL_DISABLED"), reason="Python is not free-threaded"
)
def test_importing_keeps_the_gil_disabled() -> None:
    assert not getattr(sys, "_is_gil_enabled")()  # noqa: B009
