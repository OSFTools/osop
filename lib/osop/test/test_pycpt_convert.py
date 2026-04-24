# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Tests for pycpt convert."""

import datetime

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import osop.pycpt_convert as convert


def test_standardize_dims(dataset, ds_lagged):
    """Test standardize dims with proxy array."""
    # Set up false array
    times = pd.date_range("2020-01-01", periods=2)
    time = np.array(["2020-01-01"], dtype="datetime64[ns]")
    step = pd.to_timedelta(["61 days", "92 days", "123 days"])
    forecast_month = [1, 2, 3]
    lat = [50.0, 51.0]
    lon = [-1.0, 0.0]

    out_lagged = convert._standardize_dims(
        ds_lagged,
        is_lagged=True,
        st_dim_name="start_date",
    )

    out = convert._standardize_dims(
        dataset,
        is_lagged=False,
        st_dim_name="time",
    )

    assert set(out.dims) & set(out_lagged.dims) == {"time", "step", "Y", "X"}


@pytest.fixture
def dataset():
    """Create a dataset for testing."""
    times = pd.date_range("2020-01-01", periods=2)
    time = np.array(["2020-01-01"], dtype="datetime64[ns]")
    step = pd.to_timedelta(["61 days", "92 days", "123 days"])
    forecast_month = [1, 2, 3]
    lat = [50.0, 51.0]
    lon = [-1.0, 0.0]

    dataset = xr.Dataset(
        {"var": (("time", "step", "latitude", "longitude"), np.zeros((1, 3, 2, 2)))},
        coords={
            "time": [pd.Timestamp("2020-01-01")],
            "step": pd.to_timedelta(["61 days", "92 days", "123 days"]),
            "latitude": lat,
            "longitude": lon,
        },
    )

    return dataset


@pytest.fixture
def ds_lagged():
    """Create a dataset for testing."""
    times = pd.date_range("2020-01-01", periods=2)
    time = np.array(["2020-01-01"], dtype="datetime64[ns]")
    step = pd.to_timedelta(["61 days", "92 days", "123 days"])
    forecast_month = [1, 2, 3]
    lat = [50.0, 51.0]
    lon = [-1.0, 0.0]

    ds_lagged = xr.Dataset(
        {
            "var": (
                ("start_date", "forecastMonth", "latitude", "longitude"),
                np.zeros((2, 3, 2, 2)),
            )
        },
        coords={
            "start_date": times,
            "forecastMonth": forecast_month,
            "latitude": lat,
            "longitude": lon,
        },
    )

    return ds_lagged
