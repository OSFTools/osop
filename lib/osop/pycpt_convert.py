# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Functions used to convert the ECMWF grib files to a pycpt compatible nc file.

Notes
-----
Grib formatt comes in with several restructural needs to be used for pycpt.
This set - checks that the grib exists, calculates seconds in a season frame to then convert tprate to total precip,
deals with meta handles i.e. lon/lat needs to be renamed to X/Y. Shifts lead times by 1 to align with pycpt/iri conventions
(ECMWF uses 1 to represent a iri leadtime of 0) and finally needs to shift times to sit in the middle of a month frame rather than the start.

The end result it spits out is a nc file that can be directly used by pycpt for calibrated forecasting.

"""

import logging

logger = logging.getLogger(__name__)

import datetime as dt
import errno
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from osop.util import get_tindex


def grib_open_fix(grib_file, cfgrib_kwargs=None):
    """Check grib file exists and open, with automatic lagged ensemble detection.

    For lagged ensemble datasets (UKMO), automatically reindexes to align
    multiple initialization dates into a single "start_date" dimension.

    Parameters
    ----------
    grib_file : str
        Path to grib file.
    cfgrib_kwargs : dict, optional
        Additional keyword arguments for cfgrib backend.

    Returns
    -------
    xarray.Dataset
        Opened dataset, with consistent dimension naming (time/lat/lon).
    """
    file = Path(grib_file)
    if not file.is_file():
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), grib_file)

    # Detect if this is a lagged ensemble (has indexingDate key in GRIB metadata)
    st_dim_name = get_tindex(grib_file)
    is_lagged = st_dim_name == "indexing_time"

    # Build backend kwargs and open dataset
    bkw = _build_backend_kwargs(cfgrib_kwargs, is_lagged, st_dim_name)
    ds = xr.open_dataset(grib_file, engine="cfgrib", backend_kwargs=bkw)

    # Standardize dimension names
    ds = _standardize_dims(ds, is_lagged, st_dim_name)

    # Reconstruct valid_time for lagged ensembles
    if is_lagged:
        logger.debug(f"{grib_file} detected as lagged")

        ds = _reconstruct_valid_time(ds)
    logger.info(f"{grib_file} successfully found and opened")
    return ds


def _build_backend_kwargs(cfgrib_kwargs, is_lagged, st_dim_name):
    """Build cfgrib backend kwargs for GRIB file opening.

    Parameters
    ----------
    cfgrib_kwargs : dict, optional
        User-provided keyword arguments.
    is_lagged : bool
        Whether this is a lagged ensemble.
    st_dim_name : str
        Name of the start date dimension ('time' or 'indexing_time').

    Returns
    -------
    dict
        Backend kwargs for xr.open_dataset.
    """
    bkw = {"indexpath": ""}  # Prevent creation of idx files
    if cfgrib_kwargs:
        bkw.update(cfgrib_kwargs)
    if is_lagged:
        bkw["time_dims"] = ("forecastMonth", st_dim_name)
    return bkw


def _standardize_dims(ds, is_lagged, st_dim_name):
    """Standardize dimension names for downstream processing.

    For lagged: indexing_time to time, forecastMonth to step, lat/lon to Y/X
    For burst: lat/lon to Y/X only

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to rename.
    is_lagged : bool
        Whether this is a lagged ensemble.
    st_dim_name : str
        Name of the start date dimension.

    Returns
    -------
    xarray.Dataset
        Dataset with standardized dimension names.
    """
    rename_dict = {}
    if is_lagged:
        # Lagged: indexing_time to time, forecastMonth to step
        rename_dict[st_dim_name] = "time"
        if "forecastMonth" in ds.dims:
            rename_dict["forecastMonth"] = "step"
        if "latitude" in ds.dims:
            rename_dict["latitude"] = "Y"
        if "longitude" in ds.dims:
            rename_dict["longitude"] = "X"
    else:
        # Burst: keep time as-is, just standardize lat/lon
        if "latitude" in ds.dims:
            rename_dict["latitude"] = "Y"
        if "longitude" in ds.dims:
            rename_dict["longitude"] = "X"

    if rename_dict:
        ds = ds.rename(rename_dict)
    return ds


def _reconstruct_valid_time(ds):
    """Reconstruct valid_time from time and step coordinates.

    Handles both ECMWF style (timedelta steps) and UKMO style (lead month indices).

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with time and step coordinates.

    Returns
    -------
    xarray.Dataset
        Dataset with valid_time coordinate added.
    """
    # Check if we can reconstruct valid_time: need step and a time coordinate/dim
    has_step = "step" in ds.dims and "step" in ds.coords
    has_time = "time" in ds.dims or "time" in ds.coords

    if not (has_step and has_time):
        return ds

    # Get time values (handle both dim and scalar coord cases)
    if "time" in ds.dims:
        time_vals = ds.time.values
        time_is_dim = True
    else:
        # time is a scalar coordinate
        time_vals = np.atleast_1d(ds.time.values)
        time_is_dim = False

    # Check if step is timedelta (ECMWF-like) or numeric (UKMO lead months)
    step_vals = ds.step.values

    if np.issubdtype(step_vals.dtype, np.timedelta64):
        # ECMWF style: step is timedeltas, just add directly
        if time_is_dim:
            time_2d = time_vals[:, np.newaxis]  # (time, 1)
            step_2d = step_vals[np.newaxis, :]  # (1, step)
            valid_time_2d = time_2d + step_2d
            valid_time_dims = ("time", "step")
        else:
            # For scalar time, just add the single time value to all steps
            valid_time_2d = time_vals[0] + step_vals
            valid_time_dims = ("step",)
    else:
        # UKMO style: step contains lead month indices [2, 3, 4]
        # Convert to month offsets and add to time
        if time_is_dim:
            time_2d = []
            for t in pd.to_datetime(time_vals):
                row = []
                for lead_month in step_vals:
                    # Add lead_month months to init time
                    target = t + pd.DateOffset(months=int(lead_month))
                    row.append(np.datetime64(target, "ns"))
                time_2d.append(row)
            valid_time_2d = np.array(time_2d, dtype="datetime64[ns]")
            valid_time_dims = ("time", "step")
        else:
            # For scalar time, create 1D array of valid times across steps
            row = []
            for lead_month in step_vals:
                target = pd.to_datetime(time_vals[0]) + pd.DateOffset(
                    months=int(lead_month)
                )
                row.append(np.datetime64(target, "ns"))
            valid_time_2d = np.array(row, dtype="datetime64[ns]")
            valid_time_dims = ("step",)

    # Create DataArray with proper coordinates
    ds = ds.assign_coords(valid_time=(valid_time_dims, valid_time_2d))

    return ds
