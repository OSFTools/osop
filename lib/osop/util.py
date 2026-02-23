# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""A module with utility functions for seasonal forecasts."""

import logging

import eccodes

logger = logging.getLogger(__name__)
import pandas as pd
import xarray as xr


def sel_season_time(dataset, start_year, end_year):
    """Extract data for specified years from Dec of first year to Nov of last year.

    Parameters
    ----------
    dataset : xarray.Dataset
        The data to extract the time from.
    start_year : int
        The start year of the period to extract.
    end_year : int
        The end year of the period to extract.

    Returns
    -------
    xarray.Dataset
        The data with the time extracted.
    """
    return dataset.sel(time=slice(str(start_year) + "-12-01", str(end_year) + "-11-30"))


def season_stats(dataset, start_year, end_year, stats=["mean"]):
    """Calculate the climatology mean and quantiles for a dataset.

    This is done over a given time period for standard meteorological seasons
    and the categories of the terciles for each season.

    Parameters
    ----------
    dataset : xarray.Dataset
        The data to calculate the stats on.
    start_year : int
        The start year of the period to calculate the stats.
    end_year : int
        The end year of the period to calculate the stats.
    stats : list, optional
        The statistics to calculate, default is ['mean'].
        Options: ['mean', 'seas_ts'], where seas_ts is the seasonal time series
        for all years.

    Returns
    -------
    dict
        A dictionary with the calculated stats as xarray.Datasets.
    """
    # select the time period
    dataset = sel_season_time(dataset, start_year, end_year)

    # seasonal aggregation - meaning
    # if want other aggregations will need more thought, probably using Iris
    month_length = dataset.time.dt.days_in_month

    seasonal_ds = (dataset * month_length).resample(
        time="QS-DEC"
    ).sum() / month_length.resample(time="QS-DEC").sum()

    # make an empty dict for putting the stats in
    ds_stats = {}

    for stat in stats:
        if stat == "mean":
            # calculate the seasonal mean, this is a dataset
            means = seasonal_ds.groupby("time.season").mean("time")
            ds_stats["mean"] = means
        elif stat == "seas_ts":
            ds_stats["seas_ts"] = seasonal_ds
        elif stat == "terciles":
            ds_stats["terciles"] = seasonal_ds.groupby("time.season").quantile(
                [1.0 / 3.0, 2.0 / 3.0]
            )
        else:
            raise ValueError(f"Unknown stat {stat}")

    return ds_stats


def get_tindex(infile):
    """Use eccodes to check if there is an indexing time dimension.

    Parameters
    ----------
    infile : str
        Name of file to check.

    Returns
    -------
    str
        Name of time dimension to use for indexing. Returns 'time' for burst
        ensemble and 'indexing_time' for lagged.
    """
    f = open(infile, "rb")
    gid = eccodes.codes_grib_new_from_file(f)
    key = "indexingDate"
    try:
        eccodes.codes_get(gid, key)
        st_dim_name = "indexing_time"

    except eccodes.KeyValueNotFoundError:
        st_dim_name = "time"

    eccodes.codes_release(gid)
    return st_dim_name


# move
def index(forecast_local, st_dim_name):
    """Reindex and restyle the forecast grib.

    This ensures that the data layout is consistent and compatible with
    hindcast terciles.

    Parameters
    ----------
    forecast_local : str
        File location for the grib file.
    st_dim_name : str
        Name of the start date dimension (important for lagged models).

    Returns
    -------
    xarray.Dataset
        A re-indexed x-array for forecast data.
    """
    logger.debug("Reading Forecast data from file")
    forecast_data = xr.open_dataset(
        forecast_local,
        engine="cfgrib",
        backend_kwargs=dict(time_dims=("forecastMonth", st_dim_name)),
    )
    # force dask.array using chunks on leadtime, latitude and longitude coordinate
    forecast_data = forecast_data.chunk(
        {"forecastMonth": 1, "latitude": "auto", "longitude": "auto"}
    )
    forecast_data = forecast_data.rename(
        {"latitude": "lat", "longitude": "lon", st_dim_name: "start_date"}
    )
    return forecast_data
