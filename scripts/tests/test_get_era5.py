# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Unit tests for scripts/get_era5.py."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from types import ModuleType
from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture(scope="module")
def get_era5_module():
    """Load scripts/get_era5.py as a module for testing."""
    script_path = Path(__file__).resolve().parents[1] / "get_era5.py"
    spec = spec_from_file_location("get_era5", script_path)

    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {script_path}")

    module = module_from_spec(spec)
    osop_module = ModuleType("osop")
    osop_module.__path__ = []
    pycpt_module = ModuleType("osop.pycpt_convert")
    setattr(pycpt_module, "process_grib_to_pycpt", MagicMock())
    with patch.dict(
        sys.modules,
        {"osop": osop_module, "osop.pycpt_convert": pycpt_module},
    ):
        spec.loader.exec_module(module)
    return module


def _config():
    return {
        "var": "2m_temperature",
        "hcstarty": 1993,
        "hcendy": 1994,
        "start_month": 11,
        "leads_obs": [0, 1],
        "area": [10.0, -20.0, -10.0, 20.0],
    }


def test_get_obs_skips_existing_nonempty_file(get_era5_module, tmp_path, monkeypatch):
    """Do not retrieve data when a non-empty destination already exists."""
    obs_fname = str(tmp_path / "observations.grib")
    Path(obs_fname).write_bytes(b"existing data")
    client_class = MagicMock()
    monkeypatch.setattr(get_era5_module.cdsapi, "Client", client_class)

    result = get_era5_module.get_obs(obs_fname, _config())

    assert result == obs_fname
    client_class.assert_not_called()


@pytest.mark.parametrize("existing_contents", [None, b""])
def test_get_obs_downloads_when_destination_missing_or_empty(
    get_era5_module, tmp_path, monkeypatch, existing_contents
):
    """Retrieve to a part file and rename it for missing or empty destinations."""
    obs_fname = str(tmp_path / "observations.grib")
    if existing_contents is not None:
        Path(obs_fname).write_bytes(existing_contents)

    client = MagicMock()

    def write_download(_dataset, _request, target):
        Path(target).write_bytes(b"downloaded grib data")

    client.retrieve.side_effect = write_download
    client_class = MagicMock(return_value=client)
    monkeypatch.setattr(get_era5_module.cdsapi, "Client", client_class)

    result = get_era5_module.get_obs(obs_fname, _config())

    assert result == obs_fname
    assert Path(obs_fname).read_bytes() == b"downloaded grib data"
    assert not Path(obs_fname + ".part").exists()
    client.retrieve.assert_called_once()
    request = client.retrieve.call_args.args[1]
    assert request["year"] == ["1993", "1994", "1995"]
    assert request["month"] == ["11", "12"]
    assert request["variable"] == ["2m_temperature"]


def test_get_obs_does_not_rename_when_retrieval_fails(
    get_era5_module, tmp_path, monkeypatch
):
    """Propagate retrieval errors without attempting to rename the part file."""
    obs_fname = str(tmp_path / "observations.grib")
    client = MagicMock()
    client.retrieve.side_effect = RuntimeError("download failed")
    monkeypatch.setattr(
        get_era5_module.cdsapi, "Client", MagicMock(return_value=client)
    )
    rename = MagicMock()
    monkeypatch.setattr(get_era5_module.os, "rename", rename)

    with pytest.raises(RuntimeError, match="download failed"):
        get_era5_module.get_obs(obs_fname, _config())

    rename.assert_not_called()
