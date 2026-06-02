# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Tests for regridders module."""

import os

import numpy as np
import pytest
import xarray as xr

from osop.regridders import interp_target, regrid_data_std

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")


def test_interp_target():
    """Test interp_target function."""
    domain = {"x0": -1.0, "x1": 2.0, "y0": 35.0, "y1": 39.0}
    res = 1.0
    result = interp_target(domain, res)

    # Check if the result is an xarray Dataset
    assert isinstance(result, xr.Dataset)

    # Check if the result has the expected dimensions
    assert "latitude" in result.sizes
    assert "longitude" in result.sizes

    # Check if the result has the expected coordinates
    assert np.allclose(result.latitude.values, np.arange(34.5, 39.5, res))
    assert np.allclose(result.longitude.values, np.arange(-1.5, 2.5, res))

    # Check if the result has the expected attributes
    assert "units" in result.latitude.attrs
    assert "units" in result.longitude.attrs
    assert result.latitude.attrs["units"] == "degrees_north"
    assert result.longitude.attrs["units"] == "degrees_east"


def test_interp_target_raises_x():
    """Test interp_target function raises an error for x0 >= x1."""
    domain = {"x0": 2.0, "x1": -1.0, "y0": 35.0, "y1": 39.0}
    res = 1.0
    with pytest.raises(
        ValueError, match="x0 must be smaller than x1 but got 2.0, -1.0"
    ):
        interp_target(domain, res)


def test_interp_target_raises_y():
    """Test interp_target function raises an error for y0 >= y1."""
    domain = {"x0": -1.0, "x1": 2.0, "y0": 39.0, "y1": 35.0}
    res = 1.0
    with pytest.raises(
        ValueError, match="y0 must be smaller than y1 but got 39.0, 35.0"
    ):
        interp_target(domain, res)


def test_regrid_std(xr_3x3_ds):
    """Test regrid_data_std function."""
    domain = {"x0": 1.0, "x1": 2.0, "y0": 1.0, "y1": 2.0}
    res = 1.0
    target = interp_target(domain, res)
    result = regrid_data_std(xr_3x3_ds, target)
    assert isinstance(result[0], xr.Dataset)
    assert result[1].equals(target)
    assert "latitude" in result[0].sizes
    assert "longitude" in result[0].sizes
    assert np.allclose(result[0].latitude.values, target.latitude.values)
    assert np.allclose(result[0].longitude.values, target.longitude.values)
    assert np.allclose(
        result[0].data.values, np.array([[3.0, 4.0], [6.0, 7.0]]), rtol=1e-04
    )


def test_regrid_std_fail(xr_3x3_ds):
    """Test regrid_data_std function."""
    # use an invalid target (not an xarray Dataset) to check it raises an error
    target = np.array([[1.0, 2.0], [3.0, 4.0]])
    with pytest.raises(KeyError, match="Alignment failed: please check dataset entry"):
        result = regrid_data_std(xr_3x3_ds, target)


@pytest.fixture
def xr_3x3_ds():
    """Xarray Dataset fixture with a 3x3 array.

    The data increases from top to bottom and left to right:
    [[1, 2, 3],
     [4, 5, 6],
     [7, 8, 9]]

    Coordinates:
    - latitude: top-to-bottom [0.0, 1.0, 2.0] (degrees_north)
    - longitude: left-to-right [0.0, 1.0, 2.0] (degrees_east)
    """
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).astype(float)
    lat = np.array([0.0, 1.0, 2.0])
    lon = np.array([0.0, 1.0, 2.0])

    xr_3x3_ds = xr.Dataset(
        {"data": (("latitude", "longitude"), data)},
        coords={
            "latitude": ("latitude", lat, {"units": "degrees_north"}),
            "longitude": ("longitude", lon, {"units": "degrees_east"}),
        },
    )

    return xr_3x3_ds
