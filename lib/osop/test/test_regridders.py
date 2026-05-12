# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Tests for regridders module."""

import os

import numpy as np
import pytest
import xarray as xr

from osop.regridders import interp_target, regrid_cons_masked

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")


def test_interp_target():
    """Test interp_target function."""
    domain = {"x0": -1.0, "x1": 2.0, "y0": 35.0, "y1": 39.0}
    res = 1.0
    result = interp_target(domain, res)

    # Check if the result is an xarray Dataset
    assert isinstance(result, xr.Dataset)

    # Check if the result has the expected dimensions
    assert "lat" in result.sizes
    assert "lon" in result.sizes

    # Check if the result has the expected coordinates
    assert np.allclose(result.lat.values, np.arange(34.5, 39.5, res))
    assert np.allclose(result.lon.values, np.arange(-1.5, 2.5, res))

    # Check if the result has the expected attributes
    assert "units" in result.lat.attrs
    assert "units" in result.lon.attrs
    assert result.lat.attrs["units"] == "degrees_north"
    assert result.lon.attrs["units"] == "degrees_east"


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
