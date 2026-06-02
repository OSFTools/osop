# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.


"""Functions to help with regridding for xarray."""

import logging

logger = logging.getLogger(__name__)

import numpy as np
import xarray as xr
import xarray_regrid


def interp_target(domain, res):
    """Create an interpolation target for a specific domain.

    Parameters
    ----------
    domain : dict
        Dictionary containing x0, x1, y0, y1.
    res : float
        Resolution of the target grid in degrees.

    Returns
    -------
    xarray.Dataset
        Dataset to use as interpolation target.

    Raises
    ------
    ValueError
        If x0 >= x1 or y0 >= y1.
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
            "latitude": (
                ["latitude"],
                np.arange(y0 - res / 2.0, y1 + res / 2.0, res),
                {"units": "degrees_north"},
            ),
            "longitude": (
                ["longitude"],
                np.arange(x0 - res / 2.0, x1 + res / 2.0, res),
                {"units": "degrees_east"},
            ),
        }
    )
    return target


def regrid_data_std(input_ds, target_ds):
    """Regrid dataset to match target grid resolution.

    Regrids the dataset appropriately for its type (planned expansion
    for precipitation).

    Parameters
    ----------
    input_ds : xarray.Dataset
        Data to be re-gridded.
    target_ds : xarray.Dataset
        Data set with the target grid.

    Returns
    -------
    output_ds : xarray.Dataset
        Regridded dataset to be used for analysis.
    target_ds : xarray.Dataset
        Matching target dataset (no changes).

    Raises
    ------
    KeyError
        If alignment fails due to incompatible datasets.
    """
    xr.set_options(keep_attrs=True)
    try:
        output_ds = input_ds.regrid.linear(target_ds)
    except Exception as e:
        logger.error(f"Alignment failed {e}: {e}")
        raise KeyError("Alignment failed: please check dataset entry")

    return output_ds, target_ds
