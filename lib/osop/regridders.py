"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.

Code used to help with regridding for xarray
"""

import xesmf as xe
import xarray as xr
import numpy as np


def regrid_cons_masked(source_in, var, target_in, thresh=0.5):
    """
    conservatively regrid a source data with a mask, with a tolerance
    for missing data of thresh

    Args:
                source(dataset): dataset to be regridded
        var(str):        variable name to use to find the mask
        target(dataset): dataset to use as the interpolation target
        theshold(float): threshold to use for masking, optional, default 0.5

    Returns:
                output(dataset): interpolated dataset
    """

    # make a source land-sea mask NAN = sea = 0
    # want to be able to regrid multiple times so have isel time=0
    # don't want to modify source or target so take copies
    source = source_in.copy()
    target = target_in.copy()

    source["lsm"] = xr.where(~np.isnan(source[var].isel(time=0)), 1.0, 0.0)

    regridder = xe.Regridder(source, target, "conservative_normed")

    # now regrid LSM to get a fractionsl LSM on the target grid
    lsm_out = regridder(source["lsm"], keep_attrs=True)

    # mask the source using same approach as above but integer
    source["mask"] = xr.where(~np.isnan(source[var].isel(time=0)), 1, 0)

    # for the target grid, based on where the regridded mask is at thresh or more
    target["mask"] = xr.where(lsm_out > thresh, 1.0, 0.0)

    # now regrid source respecting the new masks
    regridder = xe.Regridder(source, target, "conservative_normed")
    output = regridder(source, keep_attrs=True)

    return output


def interp_target(domain, res):
    """
    create an interpolation target for a specific domain

    Args:
        domain(dict): dictionary containing x0,x1,y0,y1
        res(float):    resolution of the target grid

    Returns:
        target(dataset): dataset to use an interpolation target
    """

    x0 = domain["x0"]
    x1 = domain["x1"]
    y0 = domain["y0"]
    y1 = domain["y1"]

    # check order
    if y1 <= y0:
        raise (ValueError(f"y0 must be smaller than y1 but got {y0}, {y1}"))
    if x1 <= x0:
        raise (ValueError(f"x0 must be smaller than x1 but got {x0}, {x1}"))

    # arange needs to be 1 unit bigger than final point, want points to be centres, so offset
    # by half a grid box
    target = xr.Dataset(
        {
            "lat": (
                ["lat"],
                np.arange(y0 - res / 2.0, y1 + res / 2.0, res),
                {"units": "degrees_north"},
            ),
            "lon": (
                ["lon"],
                np.arange(x0 - res / 2.0, x1 + res / 2.0, res),
                {"units": "degrees_east"},
            ),
        }
    )
    return target
