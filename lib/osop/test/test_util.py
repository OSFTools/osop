"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.



Tests for the util module
"""

import numpy as np
import xarray as xr
import pandas as pd
import osop.util as util
import pytest


def test_sel_season_time(dataset: xr.Dataset):
    """
    Test the sel_season_time function
    """

    start_year = 2000
    end_year = 2002

    dataset2 = util.sel_season_time(dataset, start_year, end_year)

    assert dataset2.time[0].values == np.datetime64("2000-12-01")
    assert dataset2.time[-1].values == np.datetime64("2002-11-30")


def test_season_stats_mean(dataset: xr.Dataset):
    """
    Test the season_stats function for meaning
    """

    seasonal_mean = util.season_stats(dataset, 2000, 2002)

    assert seasonal_mean["mean"].data.shape == (4, 10, 10)
    expected = [56103.76940133, 74532.95535082, 65300.0, 83650.0]
    assert np.allclose(seasonal_mean["mean"].data[:, 0, 0].data, expected)


def test_season_stats_ts(dataset: xr.Dataset):
    """
    Test the season_stats function for seasonal time-series
    """

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
    """
    Test the season_stats function for seasonal terciles
    """

    terciles = util.season_stats(dataset, 2000, 2002, stats=["terciles"])

    assert terciles["terciles"].data.shape == (4, 2, 10, 10)
    expected = [50020.436068, 62187.10273466]

    assert np.allclose(terciles["terciles"].data[0, :, 0, 0].data, expected)

    dims = [key for key in terciles["terciles"].sizes]
    assert dims == ["season", "quantile", "lat", "lon"]


def test_season_stats_foo(dataset: xr.Dataset):
    """
    Test the season_stats function for fails if given foo
    as a statistic
    """

    with pytest.raises(ValueError) as exc_info:
        util.season_stats(dataset, 2000, 2002, stats=["foo"])

    assert str(exc_info.value) == "Unknown stat foo"


@pytest.fixture
def dataset():
    """
    Create a dataset for testing
    """

    # 3 years of full calendar data with 10x10 grid
    data = np.arange(365 * 3 * 10 * 10).reshape((365 * 3, 10, 10))
    time = pd.date_range("2000-01-01", periods=365 * 3)
    dataset = xr.Dataset(
        {"data": (["time", "lat", "lon"], data)},
        coords={"time": time, "lat": np.arange(10), "lon": np.arange(10)},
    )

    return dataset
