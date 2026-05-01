# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Functions used to convert the ECMWF grib files to a pycpt compatible nc file.

Notes
-----
Currently, to use a pycpt workflow the inputting files must be standardised to a specific
and rigid format.  Note that more details on how metadata is handled by pycpt can be found here:
https://cpthelp.iri.columbia.edu/CPT_use_input_tags.html.

The overall code body of OSOP is mostly designed to work with C3S data, due to its open source
and reliable nature. However, the C3S grib files are laid out very differently from the file needed
for pycpt despite containing the same data. C3S uses these dimensions/coordinates;
number (the ensemble members),
time (initialisation date),
step (gap from int. date to forecasts in days),
surface (a height measurement),
latitude,
longitude and
valid-time (the start months of forecasts - mainly for
lagged ensembles).

Pycpt needs S (initialisation date), Ti (initial date of forecast), Tf (final date of forecast),
T (midpoint of forecast), X (latitude), Y (longitude) as well as in the case for precipitation,
accumulated rainfall instead of the tprate C3S uses. Some of these conversions are quite straight
forward - i.e. renaming latitude and longitude to X and Y. Some of the conversion is more
complicated - i.e. giving the time axis 3 coordinates and then assigning the start date,
midpoint and end date from the step coordinate in C3S. The C3S data must also have the number
coordinate averaged out of it as pycpt is not designed to work with ensemble data.
Surface must be dropped. And tprate must be converted into total accumulate precip -
 which can be done between the steps and forecast times. Is hifts lead times by 1 to align
 with pycpt/iri conventions (ECMWF uses 1 to represent a iri leadtime of 0) and finally needs
to shift times to sit in the middle of a month frame rather than the start.

The final output is a netcdf file that can be directly used by pycpt for calibrated forecasting.

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


def _accumulate_monthly(F, var_name, Sec, cast, scale=1000.0):
    """Convert precipitation to monthly totals in mm.

    F (xarray): The input dataset
    Var_name (str): The variable (precip,t2m)
    Sec (int): seconds from dataset
    cast(str): Hindcast/forecast/obs
    scale(int): meters to mm conversion

    Forecast / Hindcast:
      tprate (m s-1):  multiply by calendar-month seconds (Sec)

    ERA5 obs:
      tp (m/day): multiply by days

    Temperature:
      unchanged
    """
    if cast == "obs" and var_name == "tp":
        days = Sec / 86400.0
        monthly_mm = F * days * scale
        monthly_mm.attrs["units"] = "mm"
        return monthly_mm

    if var_name == "tprate":
        if Sec.ndim == 1:
            # Forecast: calendar-month seconds (1D)
            Sec_da = xr.DataArray(Sec, dims=("step",))
        else:
            # Hindcast: calendar-month seconds (2D)
            Sec_da = xr.DataArray(Sec, dims=("time", "step"))

        monthly_mm = F * Sec_da * scale
        monthly_mm.attrs["units"] = "mm"

        return monthly_mm

    # Temperature
    out = F
    if "units" not in out.attrs:
        out.attrs["units"] = "K"
    return out


def _aggregate_by_variable(F_sel, var_name, Sec_subset, dim="step"):
    """Aggregate forecast/hindcast data by variable: sum for precip, weighted mean for temp."""
    if var_name in ("tprate", "tp"):
        return F_sel.sum(dim, keep_attrs=True)
    elif var_name == "t2m":
        w = xr.DataArray(np.asarray(Sec_subset, dtype="float64"), dims=(dim,))
        return (F_sel * w).sum(dim, keep_attrs=True) / w.sum(dim)
    else:
        raise ValueError(f"Unknown variable: {var_name}")


def _assemble_season_coords(agg_data, T_mid, Ti_season, Tf_season, S_val):
    """Assemble aggregated data with season coordinates."""
    return agg_data.expand_dims(T=[T_mid]).assign_coords(
        Ti=("T", [Ti_season]),
        Tf=("T", [Tf_season]),
        S=("T", [S_val]),
    )


def _compute_forecast_season_bounds(Ti_val, steps_to_sum):
    """Compute season bounds from initial time for forecast."""
    season_start = pd.Timestamp(Ti_val.year, Ti_val.month, 1) - pd.DateOffset(months=1)
    Ti_season = pd.Timestamp(season_start.year, season_start.month, 1).to_datetime64()
    Tf_season = (
        season_start + pd.DateOffset(months=steps_to_sum) - pd.Timedelta(days=1)
    ).to_datetime64()
    T_mid = (Ti_season + (Tf_season - Ti_season) / np.int64(2)).astype("datetime64[ns]")
    return Ti_season, Tf_season, T_mid


def _compute_season_bounds(Ti_vals, Te_vals, Tm_vals=None):
    """Extract season boundaries and compute midpoint.

    Ti_vals (numpy.ndarray): Initial times (month starts)
    Te_vals (numpy.ndarray): End times (next month starts)
    Tm_vals (numpy.ndarray, optional): Monthly midpoint times for averaging

    Returns
    -------
    Ti_season (numpy.datetime64): Season start date
    Te_season (numpy.datetime64): Season end date (exclusive)
    Tf_season (numpy.datetime64): Season end date (inclusive)
    T_season (numpy.datetime64): Season midpoint

    """
    Ti_season = Ti_vals[0] if isinstance(Ti_vals, np.ndarray) else Ti_vals
    Te_season = Te_vals[-1] if isinstance(Te_vals, np.ndarray) else Te_vals
    Tf_season = Te_season - np.timedelta64(1, "ns")

    if Tm_vals is None:
        # Simple midpoint (forecast)
        T_season = Ti_season + (Tf_season - Ti_season) / np.int64(2)
    else:
        # Mean of monthly midpoints (hindcast)
        ints = Tm_vals.astype("datetime64[ns]").astype("int64")
        nat_min = np.iinfo(np.int64).min
        valid = ints != nat_min
        if not np.any(valid):
            raise ValueError("All monthly midpoints are NaT")
        T_season = np.datetime64(int(ints[valid].sum() // valid.sum()), "ns")

    return Ti_season, Te_season, Tf_season, T_season


def _extract_init_times(ds):
    """Extract and convert initialization times from dataset."""
    if "time" in ds.coords:
        time_vals = ds["time"].values
        if isinstance(time_vals, np.ndarray):
            return pd.to_datetime(time_vals).astype("datetime64[ns]")
        else:
            return np.array([pd.to_datetime(time_vals).asm8], dtype="datetime64[ns]")
    return np.array([np.datetime64("NaT")], dtype="datetime64[ns]")


def _finalize_season(F_season, lon_wrap="-180..180", out_var="aprod"):
    """Clean coordinates, transpose to TYX, wrap longitude, and set attributes.

    F_season (xarray.DataArray): Seasonal data array to finalize
    lon_wrap (str): Longitude wrapping style ("-180..180" or "0..360")
    out_var (str): Output variable name for renaming

    Returns
    -------
    F_season (xarray.DataArray): Finalized seasonal data array in TYX format
    """
    # Keep only PyCPT-required coords
    keep = {"T", "Ti", "Tf", "S", "Y", "X"}
    drop_these = [c for c in F_season.coords if c not in keep]
    if drop_these:
        F_season = F_season.reset_coords(names=drop_these, drop=True)

    # Check required dimensions
    for d in ("T", "Y", "X"):
        if d not in F_season.dims:
            raise ValueError(f"Missing required dimension '{d}'")

    # Reorder to TYX
    F_season = F_season.transpose("T", "Y", "X")

    # Wrap longitude if needed
    if "X" in F_season.coords:
        X = F_season["X"].values.astype(float)
        if lon_wrap == "-180..180":
            X_new = ((X + 180) % 360) - 180
        elif lon_wrap == "0..360":
            X_new = X % 360
        else:
            raise ValueError("lon_wrap must be '-180..180' or '0..360'")
        F_season = F_season.assign_coords(X=X_new).sortby("X")

    # Rename and set attributes
    F_season = F_season.rename(out_var)
    F_season.attrs.update(missing=float(-999.0))

    return F_season


def _mean_ensemble(F):
    member_dim = next(
        (d for d in ("number", "realization", "ens_member") if d in F.dims),
        None,
    )
    if member_dim:
        return F.mean(member_dim, keep_attrs=True)
    return F


def _prepare_data(ds):
    """Extract variable name and data array from dataset.

    ds (xarray.Dataset): Dataset containing data variable to process

    Returns
    -------
    F (xarray.DataArray): Processed data array
    var_name (str): Name of the variable extracted
    """
    var_name = next(iter(ds.data_vars))
    F = ds[var_name]

    return F, var_name


def _target_ym_from_init(init_date, lead_months):
    """Return target (year, month) of the first verifying month = init + lead_months."""
    init = pd.to_datetime(init_date)
    tgt = init + pd.DateOffset(months=int(lead_months))
    return int(tgt.year), int(tgt.month)


def meta_handle(
    ds,
    Ti,
    Tm,
    Te,
    Sec,
    cast,
    steps_to_sum=3,
    lead_months=1,
    lon_wrap="-180..180",
    out_var="aprod",
):
    """Process forecast, hindcast, or obs data into PyCPT-compatible seasonal format.

    Parameters
    ----------
    - ds : xarray.Dataset with climate model data
    - Ti, Tm, Te : Time coordinate arrays (initial, mid, end times)
    - Sec : Seconds per month (1D for forecast, 2D for hindcast)
    - cast : 'forecast', 'hindcast', or 'obs'
    - steps_to_sum : Number of months to combine into season
    - lead_months : Months after initialization to start season (hindcast only; forecast uses first steps)
    - out_var : Output variable name (auto-detected from var_name otherwise)
    """
    F, var_name = _prepare_data(ds)
    S_vals = _extract_init_times(ds)
    is_forecast = Ti.ndim == 1
    is_hindcast = Ti.ndim == 2

    if cast == "forecast":
        if not is_forecast:
            raise ValueError("For forecast, Ti should be 1D")

        F = F.assign_coords(
            T=("step", Tm),
            Ti=("step", Ti),
            Te=("step", Te),
        )

        S_init = S_vals[0] if len(S_vals) > 0 else np.datetime64("NaT")
        Ti_index = pd.DatetimeIndex(pd.to_datetime(Ti))
        target_year, target_month = _target_ym_from_init(S_init, lead_months)

        mask = (Ti_index.year > target_year) | (
            (Ti_index.year == target_year) & (Ti_index.month >= target_month)
        )
        j0 = int(np.where(mask)[0][0])
        j1 = j0 + steps_to_sum

        if j1 > F.sizes["step"]:
            raise ValueError(
                f"Not enough forecast steps. Need {j1}, have {F.sizes['step']}"
            )

        # Select and prepare data based on variable type
        if var_name in ("tprate", "tp"):
            monthly_mm = _accumulate_monthly(F, var_name, Sec, "forecast")
            mm_sel = monthly_mm.isel(step=slice(j0, j1))
        elif var_name == "t2m":
            mm_sel = F.isel(step=slice(j0, j1))
        else:
            raise ValueError(f"Unknown variable: {var_name}")

        # Compute season bounds and aggregate
        Ti_season, Tf_season, T_mid = _compute_forecast_season_bounds(
            Ti_index[j0], steps_to_sum
        )
        agg_data = _aggregate_by_variable(mm_sel, var_name, Sec[j0:j1])
        F_season = _assemble_season_coords(
            agg_data, T_mid, Ti_season, Tf_season, S_init
        )
        F_season = _mean_ensemble(F_season)
        out_var = "aprod" if var_name in ("tprate", "tp") else "t2m"

    elif cast == "hindcast":
        if not is_hindcast:
            raise ValueError("For hindcast, Ti should be 2D (time, step)")

        time_n, step_n = Ti.shape

        if "time" not in F.dims:
            F = F.expand_dims(time=[np.datetime64("NaT")])

        order = [d for d in ("time", "step", "number", "Y", "X") if d in F.dims]
        F = F.transpose(*order).assign_coords(
            T=(("time", "step"), Tm),
            Ti=(("time", "step"), Ti),
            Te=(("time", "step"), Te),
        )

        season_list = []

        for i in range(time_n):
            S_i = S_vals[i] if i < len(S_vals) else np.datetime64("NaT")
            init_month = int(pd.Timestamp(S_i).month)
            target_month = ((init_month - 1 + lead_months) % 12) + 1

            Ti_row = pd.DatetimeIndex(F["Ti"].isel(time=i).values).month
            try:
                j0 = np.where(Ti_row == target_month)[0][0]
            except IndexError:
                raise ValueError(
                    f"Could not find target month {target_month} for hindcast year {i} "
                    f"(init month {init_month}, lead {lead_months} months). "
                    f"Available months: {Ti_row}"
                )
            j1 = j0 + steps_to_sum

            if j1 > step_n:
                raise ValueError(
                    f"Not enough hindcast steps at year {i}. "
                    f"Need {j1} (from step {j0} + {steps_to_sum}), have {step_n}"
                )

            # Select and prepare data based on variable type
            if var_name in ("tprate", "tp"):
                monthly_mm = _accumulate_monthly(F, var_name, Sec, "hindcast")
                mm_sel_i = monthly_mm.isel(time=i, step=slice(j0, j1))
            elif var_name == "t2m":
                mm_sel_i = F.isel(time=i, step=slice(j0, j1))
            else:
                raise ValueError(f"Unknown variable: {var_name}")

            # Compute season bounds
            Tm_sel = F["T"].isel(time=i, step=slice(j0, j1)).values
            Ti_season, Te_season, Tf_season, Tmid = _compute_season_bounds(
                F["Ti"].isel(time=i, step=j0).values,
                F["Te"].isel(time=i, step=j1 - 1).values,
                Tm_sel,
            )

            # Aggregate and assemble
            Sec_subset = Sec[i, j0:j1]
            agg_data = _aggregate_by_variable(mm_sel_i, var_name, Sec_subset)
            F_season_i = _assemble_season_coords(
                agg_data, Tmid, Ti_season, Tf_season, S_i
            )
            season_list.append(F_season_i)

        F_season = xr.concat(season_list, dim="T")
        F_season = _mean_ensemble(F_season)
        out_var = "t2m" if var_name == "t2m" else "aprod"

    elif cast == "obs":
        if "time" in F.dims:
            F = F.rename({"time": "T"})
        elif "valid_time" in F.dims:
            F = F.rename({"valid_time": "T"})

        F = F.transpose("T", "Y", "X")

        F = F.assign_coords(
            T=("T", Tm.astype("datetime64[ns]")),
            Ti=("T", Ti.astype("datetime64[ns]")),
            Te=("T", Te.astype("datetime64[ns]")),
        )

        Sec_da = xr.DataArray(np.asarray(Sec, dtype="float64"), dims=("T",))
        monthly = _accumulate_monthly(F, var_name, Sec_da, "obs")

        months = pd.DatetimeIndex(monthly["Ti"].values).month
        years = pd.DatetimeIndex(monthly["Ti"].values).year

        monthly = monthly.assign_coords(
            month=("T", months),
            year=("T", years),
        )

        season_data, T_vals, Ti_vals, Tf_vals = [], [], [], []

        for y in np.unique(years):
            g = monthly.where(monthly["year"] == y, drop=True)

            if g.sizes["T"] != steps_to_sum:
                continue

            Ti_season = g["Ti"].values[0]
            Te_season = g["Te"].values[-1]
            Tf_season = Te_season - np.timedelta64(1, "ns")
            Tmid = Ti_season + (Tf_season - Ti_season) / np.int64(2)

            if var_name == "t2m":
                # temperature =  MEAN

                w = xr.DataArray(
                    (g["Te"] - g["Ti"]).astype("timedelta64[s]").astype("float64"),
                    dims=("T",),
                )

                season_val = (g * w).sum("T", keep_attrs=True) / w.sum("T")
                out_var = "t2m"

            else:
                # precipitation = SUM not accumulate monthly
                season_val = g.sum("T", keep_attrs=True)
                out_var = "prcp"

            season_data.append(season_val)
            T_vals.append(Tmid)
            Ti_vals.append(Ti_season)
            Tf_vals.append(Tf_season)

        if not season_data:
            raise ValueError("No complete seasons found in obs data")

        F_season = xr.concat(season_data, dim="T")
        F_season = F_season.assign_coords(
            T=("T", np.array(T_vals, dtype="datetime64[ns]")),
            Ti=("T", np.array(Ti_vals, dtype="datetime64[ns]")),
            Tf=("T", np.array(Tf_vals, dtype="datetime64[ns]")),
        )

    else:
        raise ValueError("cast must be 'forecast', 'hindcast', or 'obs'")

    return _finalize_season(F_season, lon_wrap, out_var)


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
