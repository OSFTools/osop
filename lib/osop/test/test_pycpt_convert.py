# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Tests for pycpt convert."""

import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# pytestmark = pytest.mark.integration
import xarray as xr

import osop.pycpt_convert as convert


@pytest.fixture
def dataset():
    """Create a dataset for testing."""
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


def test_standardize_dims(dataset, ds_lagged):
    """Test standardize dims for both lagged and burst."""
    out_burst = convert._standardize_dims(dataset, is_lagged=False, st_dim_name="time")
    # Check co-ords that are meant to be there are:
    assert set(out_burst.dims) == {"time", "step", "Y", "X"}
    # Check removal of old co-ords
    assert "latitude" not in out_burst.dims
    assert "longitude" not in out_burst.dims

    out_lagged = convert._standardize_dims(
        ds_lagged, is_lagged=True, st_dim_name="start_date"
    )
    assert set(out_lagged.dims) == {"time", "step", "Y", "X"}
    assert "start_date" not in out_lagged.dims
    assert "forecastMonth" not in out_lagged.dims
    assert "latitude" not in out_lagged.dims


def test_reconstruct_valid_time_step():
    """Test reconstruct valid time."""
    time = np.array(["2000-01-01", "2000-02-01"], dtype="datetime64[ns]")
    step = pd.to_timedelta(["0 days", "30 days"]).to_numpy()

    ds = xr.Dataset(
        {"var": (("time", "step"), np.zeros((2, 2)))},
        coords={"time": time, "step": step},
    )

    out = convert._reconstruct_valid_time(ds)

    assert "valid_time" in out.coords
    assert out["valid_time"].dims == ("time", "step")
    assert out["valid_time"].values[0, 0] == np.datetime64("2000-01-01T00:00:00", "ns")


def test_reconstruct_valid_time_scalar(dataset):
    """Test time step part 2."""
    # Reuse dataset fixture, but collapse time to a scalar coord and
    # slice step so the data variable matches the new coord length
    ds = (
        dataset[["var"]]
        .isel(latitude=0, longitude=0, time=0, step=slice(0, 2), drop=True)
        .assign_coords(
            time=np.datetime64("2020-01-01", "ns"),
            step=np.array([2, 3]),
        )
    )

    out = convert._reconstruct_valid_time(ds)

    assert "valid_time" in out.coords
    assert out["valid_time"].dims == ("step",)
    assert out["valid_time"].values[0] == np.datetime64(
        pd.Timestamp("2020-03-01"), "ns"
    )


def test_reconstruct_valid_time_timedelta_step():
    """Test valid time part3."""
    step = pd.to_timedelta(["0 days", "30 days"]).to_numpy()
    ds = xr.Dataset(
        {"var": (("step",), np.zeros((2,)))},
        coords={"time": np.datetime64("2000-01-01", "ns"), "step": step},
    )

    out = convert._reconstruct_valid_time(ds)

    assert "valid_time" in out.coords
    assert out["valid_time"].dims == ("step",)
    assert out["valid_time"].values[0] == np.datetime64("2000-01-01T00:00:00", "ns")
    assert out["valid_time"].values[1] == np.datetime64("2000-01-31T00:00:00", "ns")


def test_reconstruct_valid_time_numeric_step_time_dim_loop():
    """Test valod time branch 4."""
    time = np.array(["2000-01-01", "2000-02-01"], dtype="datetime64[ns]")
    step = np.array([2, 3, 4], dtype=np.int64)  # MUST be numeric dtype

    ds = xr.Dataset(
        {"var": (("time", "step"), np.zeros((2, 3)))},
        coords={"time": time, "step": step},
    )

    out = convert._reconstruct_valid_time(ds)

    assert "valid_time" in out.coords
    assert out["valid_time"].dims == ("time", "step")

    # First time row
    vt0 = out["valid_time"].values[0]
    assert vt0[0] == np.datetime64(pd.Timestamp("2000-03-01"), "ns")
    assert vt0[1] == np.datetime64(pd.Timestamp("2000-04-01"), "ns")
    assert vt0[2] == np.datetime64(pd.Timestamp("2000-05-01"), "ns")

    # Second time row (forces second iteration of `for t in ...`)
    vt1 = out["valid_time"].values[1]
    assert vt1[0] == np.datetime64(pd.Timestamp("2000-04-01"), "ns")


def test_build_backend_kwargs_default_non_lagged():
    """Test backend kwargs busrt."""
    out = convert._build_backend_kwargs(
        cfgrib_kwargs=None, is_lagged=False, st_dim_name="time"
    )
    assert out == {"indexpath": ""}


def test_build_backend_kwargs_merges_user_kwargs():
    """Test back end kwards."""
    out = convert._build_backend_kwargs(
        cfgrib_kwargs={"filter_by_keys": {"shortName": "t2m"}},
        is_lagged=False,
        st_dim_name="time",
    )
    assert out["indexpath"] == ""
    assert out["filter_by_keys"]["shortName"] == "t2m"


def test_build_backend_kwargs_lagged_adds_time_dims():
    out = convert._build_backend_kwargs(
        cfgrib_kwargs=None, is_lagged=True, st_dim_name="indexing_time"
    )
    assert out["time_dims"] == ("forecastMonth", "indexing_time")


@pytest.fixture(scope="module")
def grib_test_file():
    grib = (
        Path(__file__).parent
        / "test_data"
        / "test_pycptf"
        / "pycpt_ecmwf_51_1993-2016_monthly_mean_5_234_40:10:10:40_total_precipitation.grib"
    )

    assert grib.is_file(), f"Test GRIB not found: {grib}"
    return grib


def test_grib_open_fix_real_grib(grib_test_file):
    ds = convert.grib_open_fix(grib_test_file)

    # PyCPT contract checks
    assert set(ds.dims) >= {"time", "step", "Y", "X"}
    assert "valid_time" in ds.coords

    # basic sanity checks
    assert ds.sizes["Y"] > 1
    assert ds.sizes["X"] > 1
    assert ds.sizes["step"] > 0
