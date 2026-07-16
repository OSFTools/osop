# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Unit tests for scripts/get_chirps.py."""

from importlib.util import module_from_spec, spec_from_file_location
import logging
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, mock_open, patch

import pytest
import requests


@pytest.fixture(scope="module")
def get_chirps_module():
    """Load scripts/get_chirps.py as a module for testing."""
    script_path = Path(__file__).resolve().parents[1] / "get_chirps.py"
    spec = spec_from_file_location("get_chirps", script_path)

    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {script_path}")

    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parse_args(get_chirps_module, monkeypatch):
    """Test parse_args returns expected values from CLI input."""
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "get_chirps.py",
            "--month",
            "3",
            "--leads",
            "1,2,3",
            "--area",
            "10,-20,-10,20",
            "--downloaddir",
            "/tmp/downloads",
            "--logdir",
            "/tmp/logs",
            "--pycptdir",
            "/tmp/pycpt",
            "--pycpt",
            "False",
            "--years",
            "1990,1995",
        ],
    )

    args = get_chirps_module.parse_args()

    assert args.month == "3"
    assert args.leads == "1,2,3"
    assert args.area == "10,-20,-10,20"
    assert args.downloaddir == "/tmp/downloads"
    assert args.logdir == "/tmp/logs"
    assert args.pycptdir == "/tmp/pycpt"
    assert args.pycpt == "False"
    assert args.years == "1990,1995"


def test_unpack_args_and_run_with_years(get_chirps_module, monkeypatch):
    """Test unpack_args_and_run branch when explicit years are provided."""
    captured = {}

    def fake_get_obs(downloaddir, config):
        captured["downloaddir"] = downloaddir
        captured["config"] = config

    monkeypatch.setattr(get_chirps_module, "get_obs", fake_get_obs)
    monkeypatch.setattr(get_chirps_module.logging, "basicConfig", lambda **kwargs: None)

    args = SimpleNamespace(
        month="3",
        leads="1,3,6",
        area="10,-20,-10,20",
        downloaddir="/tmp/downloads",
        logdir="/tmp/logs",
        pycptdir="/tmp/pycpt",
        pycpt="False",
        years="1990,1995",
    )

    get_chirps_module.unpack_args_and_run(args)

    assert captured["downloaddir"] == "/tmp/downloads"
    assert captured["config"]["month"] == 3
    assert captured["config"]["leadtime_month"] == [0, 2, 5]
    assert captured["config"]["area"] == [10.0, -20.0, -10.0, 20.0]
    assert captured["config"]["hcstarty"] == 1990
    assert captured["config"]["hcendy"] == 1995


def test_unpack_args_and_run_uses_default_years(get_chirps_module, monkeypatch):
    """Test unpack_args_and_run branch when years are not provided."""
    captured = {}

    def fake_get_obs(downloaddir, config):
        captured["downloaddir"] = downloaddir
        captured["config"] = config

    monkeypatch.setattr(get_chirps_module, "get_obs", fake_get_obs)
    monkeypatch.setattr(get_chirps_module.logging, "basicConfig", lambda **kwargs: None)

    args = SimpleNamespace(
        month="11",
        leads="2,4",
        area="5,6,7,8",
        downloaddir="/tmp/downloads",
        logdir="/tmp/logs",
        pycptdir="/tmp/pycpt",
        pycpt="False",
        years=None,
    )

    get_chirps_module.unpack_args_and_run(args)

    assert captured["downloaddir"] == "/tmp/downloads"
    assert captured["config"]["month"] == 11
    assert captured["config"]["leadtime_month"] == [1, 3]
    assert captured["config"]["area"] == [5.0, 6.0, 7.0, 8.0]
    assert captured["config"]["hcstarty"] == 1993
    assert captured["config"]["hcendy"] == 2016


def test_unpack_args_and_run_raises_for_pycpt_true(get_chirps_module, monkeypatch):
    """Test unpack_args_and_run raises NotImplementedError when pycpt is True."""
    state = {"calls": 0}

    def fake_get_obs(downloaddir, config):
        state["calls"] += 1

    monkeypatch.setattr(get_chirps_module, "get_obs", fake_get_obs)
    monkeypatch.setattr(get_chirps_module.logging, "basicConfig", lambda **kwargs: None)

    args = SimpleNamespace(
        month="7",
        leads="1,2",
        area="1,2,3,4",
        downloaddir="/tmp/downloads",
        logdir="/tmp/logs",
        pycptdir="/tmp/pycpt",
        pycpt="True",
        years=None,
    )

    with pytest.raises(
        NotImplementedError, match="pycpt calibration not yet implemented"
    ):
        get_chirps_module.unpack_args_and_run(args)

    # get_obs is called before the pycpt not-implemented guard is raised.
    assert state["calls"] == 1


def test_get_obs_downloads_and_writes_chunks(get_chirps_module):
    """Test that get_obs downloads each month and writes response chunks."""
    response = MagicMock()
    response.status_code = 200
    response.raise_for_status.return_value = None
    response.iter_content.return_value = [b"abc", b"", b"def"]

    config = {
        "hcstarty": 2020,
        "hcendy": 2020,
        "month": 1,
        "leadtime_month": [0, 1],
        "area": [10.0, -20.0, -10.0, 20.0],
    }

    with (
        patch.object(get_chirps_module.Path, "exists", return_value=False),
        patch.object(
            get_chirps_module.requests, "get", return_value=response
        ) as get_mock,
        patch("builtins.open", mock_open()) as open_mock,
    ):
        get_chirps_module.get_obs("/tmp/chirps", config)

    assert get_mock.call_count == 2
    assert open_mock.call_count == 2
    assert open_mock().write.call_count == 4

    # Verify exact URL and kwargs passed to requests.get
    expected_urls = [
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2020.01.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2020.02.tif",
    ]
    observed_urls = [call.args[0] for call in get_mock.call_args_list]
    assert observed_urls == expected_urls
    for call in get_mock.call_args_list:
        assert call.kwargs["timeout"] == 30
        assert call.kwargs["stream"] is True


def test_get_obs_skips_when_file_exists(get_chirps_module):
    """Test that get_obs does not download or write if file already exists."""
    config = {
        "hcstarty": 2020,
        "hcendy": 2020,
        "month": 1,
        "leadtime_month": [0, 1],
        "area": [10.0, -20.0, -10.0, 20.0],
    }

    with (
        patch.object(get_chirps_module.Path, "exists", return_value=True),
        patch.object(get_chirps_module.requests, "get") as get_mock,
        patch("builtins.open", mock_open()) as open_mock,
    ):
        get_chirps_module.get_obs("/tmp/chirps", config)

    get_mock.assert_not_called()
    open_mock.assert_not_called()


def test_get_obs_does_not_write_on_request_error(get_chirps_module):
    """Test that get_obs handles request errors and avoids writing files."""
    response = MagicMock()
    response.status_code = 403
    response.raise_for_status.side_effect = requests.exceptions.HTTPError(
        "403 Client Error",
        response=response,
    )

    config = {
        "hcstarty": 2020,
        "hcendy": 2020,
        "month": 1,
        "leadtime_month": [0, 1],
        "area": [10.0, -20.0, -10.0, 20.0],
    }

    with (
        patch.object(get_chirps_module.Path, "exists", return_value=False),
        patch.object(
            get_chirps_module.requests, "get", return_value=response
        ) as get_mock,
        patch("builtins.open", mock_open()) as open_mock,
        pytest.raises(
            RuntimeError, match="Finished downloading CHIRPS data with 2 failures"
        ),
    ):
        get_chirps_module.get_obs("/tmp/chirps", config)

    assert get_mock.call_count == 2
    open_mock.assert_not_called()

    # Verify exact URL and kwargs passed to requests.get
    expected_urls = [
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2020.01.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2020.02.tif",
    ]
    observed_urls = [call.args[0] for call in get_mock.call_args_list]
    assert observed_urls == expected_urls
    for call in get_mock.call_args_list:
        assert call.kwargs["timeout"] == 30
        assert call.kwargs["stream"] is True


def test_get_obs_clamps_start_year_before_1981(get_chirps_module, caplog):
    """Test that years before 1981 are clamped to 1981."""
    response = MagicMock()
    response.status_code = 200
    response.raise_for_status.return_value = None
    response.iter_content.return_value = [b"ok"]

    config = {
        "hcstarty": 1979,
        "hcendy": 1981,
        "month": 1,
        "leadtime_month": [0, 1],
        "area": [10.0, -20.0, -10.0, 20.0],
    }

    with caplog.at_level(logging.WARNING):
        with (
            patch.object(get_chirps_module.Path, "exists", return_value=False),
            patch.object(
                get_chirps_module.requests,
                "get",
                return_value=response,
            ) as get_mock,
            patch("builtins.open", mock_open()) as open_mock,
        ):
            get_chirps_module.get_obs("/tmp/chirps", config)

    assert get_mock.call_count == 2
    assert open_mock.call_count == 2
    assert (
        "Data from before 1981 not available, setting start year to 1981" in caplog.text
    )

    # Verify exact URL and kwargs passed to requests.get
    expected_urls = [
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.1981.01.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.1981.02.tif",
    ]
    observed_urls = [call.args[0] for call in get_mock.call_args_list]
    assert observed_urls == expected_urls
    for call in get_mock.call_args_list:
        assert call.kwargs["timeout"] == 30
        assert call.kwargs["stream"] is True


def test_get_obs_clamps_end_year_after_current_year(get_chirps_module, caplog):
    """Test that years after the current year are clamped to the current year."""
    response = MagicMock()
    response.status_code = 200
    response.raise_for_status.return_value = None
    response.iter_content.return_value = [b"ok"]

    config = {
        "hcstarty": 2020,
        "hcendy": 3000,
        "month": 1,
        "leadtime_month": [0, 1],
        "area": [10.0, -20.0, -10.0, 20.0],
    }

    real_datetime = get_chirps_module.datetime

    class FixedDateTime(real_datetime):
        @classmethod
        def now(cls, tz=None):
            return real_datetime(2024, 1, 1)

    with caplog.at_level(logging.WARNING):
        with (
            patch.object(get_chirps_module, "datetime", FixedDateTime),
            patch.object(get_chirps_module.Path, "exists", return_value=False),
            patch.object(
                get_chirps_module.requests,
                "get",
                return_value=response,
            ) as get_mock,
            patch("builtins.open", mock_open()) as open_mock,
        ):
            get_chirps_module.get_obs("/tmp/chirps", config)

    assert get_mock.call_count == 10
    assert open_mock.call_count == 10
    assert "Data from after 2024 not available, setting end year to 2024" in caplog.text

    # Verify exact URL and kwargs passed to requests.get
    expected_urls = [
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2020.01.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2020.02.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2021.01.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2021.02.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2022.01.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2022.02.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2023.01.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2023.02.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2024.01.tif",
        "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.2024.02.tif",
    ]
    observed_urls = [call.args[0] for call in get_mock.call_args_list]
    assert observed_urls == expected_urls
    for call in get_mock.call_args_list:
        assert call.kwargs["timeout"] == 30
        assert call.kwargs["stream"] is True
