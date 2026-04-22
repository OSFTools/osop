# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Tests for pycpt convert."""

import numpy as np
import pandas as pd
import xarray as xr
import pytest

from osop.pycpt_convert import calculate_month_metrics, choose_month_starts, _standardize_dims


def test_calculate_month_metrics_basic_month():
    """Test caululate month metrics.
    should assert mid point and duration of a calendar time frame.
    """
    month_start = pd.Timestamp("2001-01-01")
    next_start = pd.Timestamp("2001-02-01")

    mid, sec = calculate_month_metrics(month_start, next_start)

    assert mid == np.datetime64("2001-01-16T12:00:00")
    assert sec == 31 * 86400


def test_choose_month_starts():
    """Test choose month starts w/ condition 1 and 2."""
    VT_start = np.array(
        ["2001-01-01", "2001-02-01", "2001-03-01"],
        dtype="datetime64[ns]",
    )
    VT_start_next = np.array(
        ["2001-02-01", "2001-03-01", "2001-04-01"],
        dtype="datetime64[ns]",
    )
    VT_start_prev = np.array(
        ["2000-12-01", "2001-01-01", "2001-02-01"],
        dtype="datetime64[ns]",
    )

    i_slice = slice(1, 3)  # Feb, Mar

    month_start, next_start = choose_month_starts(
        i_slice, VT_start, VT_start_next, VT_start_prev, off=1
    )

    i_slice_2 = slice(0, 2)  # Jan, Feb

    month_start_2, next_start_2 = choose_month_starts(
        i_slice_2, VT_start, VT_start_next, VT_start_prev, off=2
    )

    assert np.all(month_start == VT_start[i_slice])
    assert np.all(next_start == VT_start_next[i_slice])
    assert np.all(month_start_2 == VT_start_prev[i_slice_2])
    assert np.all(next_start_2 == VT_start[i_slice_2])


def test_choose_month_starts_invalid_offset():
    """Test choose_month_starts with error force."""
    with pytest.raises(RuntimeError, match="offset must be 1 or 2"):
        choose_month_starts(
            slice(0, 1),
            np.array([], dtype="datetime64[ns]"),
            np.array([], dtype="datetime64[ns]"),
            np.array([], dtype="datetime64[ns]"),
            off=99,
        )



def test_standardize_dims(): 
    """Test standardize dims with proxy array."""
    #Set up false array
    times = pd.date_range("2020-01-01", periods=2)
    forecast_month = [1, 2, 3]
    lat = [50.0, 51.0]
    lon = [-1.0, 0.0]

    ds_lagged = xr.Dataset(
        {
            "var": (("start_date", "forecastMonth", "latitude", "longitude"),
                    np.zeros((2, 3, 2, 2)))
        },
        coords={
            "start_date": times,
            "forecastMonth": forecast_month,
            "latitude": lat,
            "longitude": lon,
        },
    )


    out_lagged = _standardize_dims(
            ds_lagged,
            is_lagged=True,
            st_dim_name="start_date",
        )
    
    out = _standardize_dims(
            ds_lagged,
            is_lagged=False,
            st_dim_name="start_date",
        )

    assert set(out.dims) == {"time", "step", "Y", "X"}



