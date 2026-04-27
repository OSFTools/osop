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


def test_mean_ensemble():
    """Test Mean ensemble - average the ensembel for each coord."""
    # Set up false arrays
    valid_time = pd.to_datetime(["2020-01-01"])
    forecast_month = [1, 2, 3]
    lat = [50.0, 51.0]
    lon = [-1.0, 0.0]
    number = [2, 3]

    ds_lagged = xr.Dataset(
        {
            "var": (
                ("valid_time", "forecastMonth", "latitude", "longitude", "number"),
                np.zeros((1, 3, 2, 2, 2)),
            )
        },
        coords={
            "valid_time": valid_time,
            "forecastMonth": forecast_month,
            "latitude": lat,
            "longitude": lon,
            "number": number,
        },
    )

    da_out = convert._mean_ensemble(ds_lagged["var"])

    np.testing.assert_allclose(da_out.values, 0.0)
    assert "number" not in da_out
