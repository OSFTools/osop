# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Unit tests for scripts/get_chirps.py."""

from importlib.util import module_from_spec, spec_from_file_location
import logging
from pathlib import Path
import sys
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


def test_get_obs_downloads_and_writes_chunks(get_chirps_module):
    """Test that get_obs downloads each month and writes response chunks."""
    response = MagicMock()
    response.status_code = 200
    response.raise_for_status.return_value = None
    response.iter_content.return_value = [b"abc", b"", b"def"]

    config = {"hcstarty": 2020, "hcendy": 2020}

    with (
        patch.object(get_chirps_module.Path, "exists", return_value=False),
        patch.object(
            get_chirps_module.requests, "get", return_value=response
        ) as get_mock,
        patch("builtins.open", mock_open()) as open_mock,
    ):
        get_chirps_module.get_obs("/tmp/chirps", config)

    assert get_mock.call_count == 12
    assert open_mock.call_count == 12
    assert open_mock().write.call_count == 24


def test_get_obs_skips_when_file_exists(get_chirps_module):
    """Test that get_obs does not download or write if file already exists."""
    config = {"hcstarty": 2020, "hcendy": 2020}

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

    config = {"hcstarty": 2020, "hcendy": 2020}

    with (
        patch.object(get_chirps_module.Path, "exists", return_value=False),
        patch.object(
            get_chirps_module.requests, "get", return_value=response
        ) as get_mock,
        patch("builtins.open", mock_open()) as open_mock,
    ):
        get_chirps_module.get_obs("/tmp/chirps", config)

    # Some implementations retry 403 responses once, so accept >= 12 calls.
    assert get_mock.call_count >= 12
    open_mock.assert_not_called()


def test_get_obs_clamps_start_year_before_1981(get_chirps_module, caplog):
    """Test that years before 1981 are clamped to 1981."""
    response = MagicMock()
    response.status_code = 200
    response.raise_for_status.return_value = None
    response.iter_content.return_value = [b"ok"]

    config = {"hcstarty": 1979, "hcendy": 1981}

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

    assert get_mock.call_count == 12
    assert open_mock.call_count == 12
    assert (
        "Data from before 1981 not available, setting start year to 1981" in caplog.text
    )


def test_get_obs_clamps_end_year_after_current_year(get_chirps_module, caplog):
    """Test that years after the current year are clamped to the current year."""
    response = MagicMock()
    response.status_code = 200
    response.raise_for_status.return_value = None
    response.iter_content.return_value = [b"ok"]

    config = {"hcstarty": 2024, "hcendy": 2026}

    mocked_now = MagicMock()
    mocked_now.year = 2024

    with caplog.at_level(logging.WARNING):
        with (
            patch.object(get_chirps_module, "datetime") as datetime_mock,
            patch.object(get_chirps_module.Path, "exists", return_value=False),
            patch.object(
                get_chirps_module.requests,
                "get",
                return_value=response,
            ) as get_mock,
            patch("builtins.open", mock_open()) as open_mock,
        ):
            datetime_mock.now.return_value = mocked_now
            get_chirps_module.get_obs("/tmp/chirps", config)

    assert get_mock.call_count == 12
    assert open_mock.call_count == 12
    assert "Data from after 2024 not available, setting end year to 2024" in caplog.text
