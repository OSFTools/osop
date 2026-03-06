# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Tests for ens_plotting module."""

import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison
import numpy as np
import pytest
import xarray as xr

from osop import ens_plotting


def test_truncate_colormap():
    """Test truncate_colormap function."""
    cmap = ens_plotting.truncate_colormap(plt.cm.Greys, minval=0.2)
    assert cmap.name == "trunc(Greys,0.20,1.00)"


def test_get_cmap_defaults():
    """Test get_cmap function with defaults (precip_cs=False, wmo_cs=True)."""
    cmap_below, cmap_normal, cmap_above = ens_plotting.get_cmap()
    names = [cmap_below.name, cmap_normal.name, cmap_above.name]
    assert names == ["trunc(Blues,0.20,1.00)", "trunc(Greys,0.20,1.00)", "YlOrRd"]


def test_get_cmap_wmo_cs_precip():
    """Test get_cmap function for precip with WMO colours."""
    cmap_below, cmap_normal, cmap_above = ens_plotting.get_cmap(
        precip_cs=True, wmo_cs=True
    )
    names = [cmap_below.name, cmap_normal.name, cmap_above.name]
    assert names == ["YlOrRd", "trunc(Greys,0.20,1.00)", "trunc(Greens,0.20,1.00)"]


def test_get_cmap_nwmo_temp():
    """Test get_cmap function for not precip and not WMO colors."""
    cmap_below, cmap_normal, cmap_above = ens_plotting.get_cmap(
        precip_cs=False, wmo_cs=False
    )
    names = [cmap_below.name, cmap_normal.name, cmap_above.name]
    assert names == [
        "trunc(Blues,0.20,1.00)",
        "trunc(Greys,0.20,1.00)",
        "trunc(Reds,0.20,1.00)",
    ]


def test_get_cmap_wmo_temp():
    """Test get_cmap function for precip and non WMO colors."""
    cmap_below, cmap_normal, cmap_above = ens_plotting.get_cmap(
        precip_cs=True, wmo_cs=False
    )
    names = [cmap_below.name, cmap_normal.name, cmap_above.name]
    assert names == [
        "trunc(BrBG,0.00,0.15)_r",
        "trunc(Greys,0.20,1.00)",
        "trunc(Purples,0.50,1.00)",
    ]


@pytest.fixture
def mme():
    # create dummy data. 3x3 grid with 3 members
    # below starts at 0, up to 0.8
    b = np.linspace(0.0, 0.8, 9).reshape(3, 3)

    # normal starts at 0.6 then down to 0.1
    n = np.append(np.linspace(0.6, 0.1, 6), np.zeros(3)).reshape(3, 3)

    # also want 1 pixel with all ~0.33
    b[1, 1] = 0.33
    n[1, 1] = 0.33

    # calculate the above category so that all sum to 1
    a = np.ones((3, 3)) - b - n

    # package into 3x3x3 array and scale to 100 percent
    precipitation = np.stack([b, n, a], axis=0) * 100.0

    # dimensions
    C = np.array([1, 2, 3])
    Y = np.linspace(15.0, 35.0, 3)
    X = np.linspace(15.0, 35.0, 3)

    mme = xr.Dataset(
        data_vars=dict(
            precipitation=(["C", "Y", "X"], precipitation),
        ),
        coords=dict(X=X, Y=Y, C=C),
    )

    return mme


@pytest.fixture
def mask():
    # create dummy mask
    data = np.ones((3, 3))
    data[2, 2] = np.nan
    Y = np.linspace(15.0, 35.0, 3)
    X = np.linspace(15.0, 35.0, 3)

    mask = xr.DataArray(data=data, dims=["Y", "X"], coords=dict(X=X, Y=Y))
    return mask


@image_comparison(
    baseline_images=["terciles"],
    tol=1.0,
    remove_text=True,
    extensions=["png"],
    style="mpl20",
)
def test_plot_tercile_fc(mme):
    """Test plot_tercile_fc function."""
    atitle = "Test plot"
    shpfile = cfeature.NaturalEarthFeature(
        category="cultural", name="admin_0_countries", scale="10m", facecolor="none"
    )
    fig = ens_plotting.plot_tercile_fc(mme, atitle, map_setting=shpfile)


@image_comparison(
    baseline_images=["terciles_mask"],
    tol=1.0,
    remove_text=True,
    extensions=["png"],
    style="mpl20",
)
def test_plot_tercile_fc_mask(mme, mask):
    """Test plot_tercile_fc function with mask, no coastlines, plotting temperature."""
    atitle = "Test plot"
    mme = mme.rename({"precipitation": "temperature"})
    fig = ens_plotting.plot_tercile_fc(
        mme, atitle, mask=mask, var="temperature", map_setting="False"
    )
