# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Unit tests for scripts/get_any_hindcast.py."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from unittest.mock import MagicMock

import pytest


@pytest.fixture(scope="module")
def get_any_hindcast_module():
    """Load scripts/get_any_hindcast.py as a module for testing."""
    script_path = Path(__file__).resolve().parents[1] / "get_any_hindcast.py"
    spec = spec_from_file_location("get_any_hindcast", script_path)

    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {script_path}")

    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _call_arguments(tmp_path, pycpt, years):
    return (
        "ecmwf",
        5,
        3,
        [1, 2],
        [10.0, -20.0, -10.0, 20.0],
        "10.0.-20.0.-10.0.20.0",
        "2m_temperature",
        str(tmp_path),
        pycpt,
        years,
    )


def _filename(tmp_path, pycpt, years):
    if years == "hc":
        years = [str(year) for year in range(1993, 2017)]
    prefix = "pycpt_" if pycpt == "download_on" else ""
    return (
        tmp_path / f"{prefix}ecmwf_5_{years[0]}-{years[-1]}_monthly_mean_3_12_"
        "10.0.-20.0.-10.0.20.0_2m_temperature.grib"
    )


@pytest.mark.parametrize("pycpt", ["download_on", "False"])
def test_do_cdsapi_call_skips_existing_nonempty_file(
    get_any_hindcast_module, tmp_path, monkeypatch, pycpt
):
    """Skip retrieval for an existing non-empty file in either mode."""
    fname = _filename(tmp_path, pycpt, "hc")
    fname.write_bytes(b"existing data")
    client_class = MagicMock()
    monkeypatch.setattr(get_any_hindcast_module.cdsapi, "Client", client_class)

    get_any_hindcast_module.do_cdsapi_call(*_call_arguments(tmp_path, pycpt, "hc"))

    client_class.assert_not_called()


@pytest.mark.parametrize("pycpt", ["download_on", "False"])
@pytest.mark.parametrize("existing_contents", [None, b""])
@pytest.mark.parametrize(
    ("years", "expected_years"),
    [
        ("hc", [str(year) for year in range(1993, 2017)]),
        ([2001, 2002], [2001, 2002]),
    ],
)
def test_do_cdsapi_call_downloads_missing_or_empty_file(
    get_any_hindcast_module,
    tmp_path,
    monkeypatch,
    pycpt,
    existing_contents,
    years,
    expected_years,
):
    """Download and rename data for a missing or zero-byte file."""
    fname = _filename(tmp_path, pycpt, years)
    if existing_contents is not None:
        fname.write_bytes(existing_contents)

    client = MagicMock()

    def write_download(_dataset, _request, target):
        Path(target).write_bytes(b"downloaded grib data")

    client.retrieve.side_effect = write_download
    monkeypatch.setattr(
        get_any_hindcast_module.cdsapi,
        "Client",
        MagicMock(return_value=client),
    )

    get_any_hindcast_module.do_cdsapi_call(*_call_arguments(tmp_path, pycpt, years))

    assert fname.read_bytes() == b"downloaded grib data"
    assert not list(tmp_path.glob(f"{fname.name}.*.part"))
    client.retrieve.assert_called_once()
    dataset, request, target = client.retrieve.call_args.args
    assert dataset == "seasonal-monthly-single-levels"
    assert request == {
        "format": "grib",
        "originating_centre": "ecmwf",
        "system": 5,
        "variable": ["2m_temperature"],
        "product_type": "monthly_mean",
        "year": expected_years,
        "month": "03",
        "leadtime_month": [1, 2],
        "area": [10.0, -20.0, -10.0, 20.0],
    }
    assert Path(target).parent == fname.parent
    assert target.endswith(".part")
    assert Path(target).name.startswith(f"{fname.name}.")
