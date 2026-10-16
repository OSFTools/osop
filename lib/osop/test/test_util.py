# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Tests for the util module."""

from concurrent.futures import ThreadPoolExecutor
from threading import Event

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import osop.util as util
from osop.util import download_and_rename_cds


def test_sel_season_time(dataset: xr.Dataset):
    """Test the sel_season_time function."""
    start_year = 2000
    end_year = 2002

    dataset2 = util.sel_season_time(dataset, start_year, end_year)

    assert dataset2.time[0].values == np.datetime64("2000-12-01")
    assert dataset2.time[-1].values == np.datetime64("2002-11-30")


def test_season_stats_mean(dataset: xr.Dataset):
    """Test the season_stats function for meaning."""
    seasonal_mean = util.season_stats(dataset, 2000, 2002)

    assert seasonal_mean["mean"].data.shape == (4, 10, 10)
    expected = [56103.76940133, 74532.95535082, 65300.0, 83650.0]
    assert np.allclose(seasonal_mean["mean"].data[:, 0, 0].data, expected)


def test_season_stats_ts(dataset: xr.Dataset):
    """Test the season_stats function for seasonal time-series."""
    seasonal_ts = util.season_stats(dataset, 2000, 2002, stats=["seas_ts"])

    assert seasonal_ts["seas_ts"].data.shape == (8, 10, 10)
    expected = [
        37853.76940133,
        47050.0,
        56282.95535082,
        65400.0,
        74353.76940133,
        83550.0,
        92782.95535082,
        101900.0,
    ]

    assert np.allclose(seasonal_ts["seas_ts"].data[:, 0, 0].data, expected)


def test_season_tercile(dataset: xr.Dataset):
    """Test the season_stats function for seasonal terciles."""
    terciles = util.season_stats(dataset, 2000, 2002, stats=["terciles"])

    assert terciles["terciles"].data.shape == (4, 2, 10, 10)
    expected = [50020.436068, 62187.10273466]

    assert np.allclose(terciles["terciles"].data[0, :, 0, 0].data, expected)

    dims = [key for key in terciles["terciles"].sizes]
    assert dims == ["season", "quantile", "lat", "lon"]


def test_season_stats_foo(dataset: xr.Dataset):
    """Test the season_stats function for fails if given foo as a statistic."""
    with pytest.raises(ValueError) as exc_info:
        util.season_stats(dataset, 2000, 2002, stats=["foo"])

    assert str(exc_info.value) == "Unknown stat foo"


@pytest.fixture
def dataset():
    """Create a dataset for testing."""
    # 3 years of full calendar data with 10x10 grid
    data = np.arange(365 * 3 * 10 * 10).reshape((365 * 3, 10, 10))
    time = pd.date_range("2000-01-01", periods=365 * 3)
    dataset = xr.Dataset(
        {"data": (["time", "lat", "lon"], data)},
        coords={"time": time, "lat": np.arange(10), "lon": np.arange(10)},
    )

    return dataset


def make_mock_client_success(tmp_path):
    class MockClient:
        def retrieve(self, dataset, request, target_path):
            with open(target_path, "wb") as f:
                f.write(b"data")

    return MockClient()


def make_mock_client_failure(tmp_path):
    class MockClient:
        def retrieve(self, dataset, request, target_path):
            with open(target_path, "wb") as f:
                f.write(b"partial")
            raise RuntimeError("cdsapi error")

    return MockClient()


def test_download_and_rename_cds_success(tmp_path):
    client = make_mock_client_success(tmp_path)
    final = tmp_path / "out.grib"

    result = download_and_rename_cds(client, "dataset", {"a": 1}, str(final))

    assert result == str(final)
    assert final.exists()
    assert not list(tmp_path.glob("out.grib.*.part"))
    with open(final, "rb") as f:
        assert f.read() == b"data"


def test_download_and_rename_cds_failure_cleans_partial(tmp_path):
    client = make_mock_client_failure(tmp_path)
    final = tmp_path / "out2.grib"

    with pytest.raises(RuntimeError):
        download_and_rename_cds(client, "dataset", {"a": 1}, str(final))

    assert not list(tmp_path.glob("out2.grib.*.part"))
    assert not final.exists()


def test_download_and_rename_cds_overlapping_downloads_use_distinct_temp_files(
    tmp_path,
):
    """Keep one failed download from publishing another call's partial data."""
    final = tmp_path / "overlapping.grib"
    first_started = Event()
    allow_first_failure = Event()

    class OverlappingClient:
        def retrieve(self, dataset, request, target_path):
            if request["download"] == "first":
                with open(target_path, "wb") as file:
                    file.write(b"partial data")
                first_started.set()
                if not allow_first_failure.wait(timeout=5):
                    raise TimeoutError("first download was not released")
                raise RuntimeError("first download failed")

            if not first_started.wait(timeout=5):
                raise TimeoutError("first download did not start")

    client = OverlappingClient()
    with ThreadPoolExecutor(max_workers=2) as executor:
        first = executor.submit(
            download_and_rename_cds,
            client,
            "dataset",
            {"download": "first"},
            str(final),
        )
        assert first_started.wait(timeout=5)
        second = executor.submit(
            download_and_rename_cds,
            client,
            "dataset",
            {"download": "second"},
            str(final),
        )
        try:
            with pytest.raises(IOError, match="produced empty file"):
                second.result(timeout=5)
        finally:
            allow_first_failure.set()

        with pytest.raises(RuntimeError, match="first download failed"):
            first.result(timeout=5)

    assert not final.exists()
    assert not list(tmp_path.glob("overlapping.grib.*.part"))
