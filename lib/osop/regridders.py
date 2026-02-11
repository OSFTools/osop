# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.


"""Functions to help with regridding for xarray."""

import numpy as np
import xarray as xr
import xesmf as xe


def regrid_cons_masked(source_in, var, target_in, thresh=0.5):
    """Conservatively regrid a source data with a mask, with a tolerance for missing data of thresh.

    Parameters
    ----------
    source_in : dataset
        Dataset to be regridded
    var : str
        Variable name to use to find the mask
    target_in : dataset
        Dataset to use as the interpolation target
    thresh : float, optional
        Threshold to use for masking, default 0.5

    Returns
    -------
    output : dataset
        Interpolated dataset
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
    """Create an interpolation target for a specific domain.

    Parameters
    ----------
    domain : dict
        Dictionary containing x0,x1,y0,y1
    res : float
        Resolution of the target grid in degrees

    Returns
    -------
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
