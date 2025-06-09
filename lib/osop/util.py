"""
A module with utility functions for seasonal FCs
"""

import xarray as xr


def sel_season_time(dataset, start_year, end_year):
    """Extract the data for the specified years starting with Dec of
    first year and ending with Nov of last year

    Parameters:
    dataset: xarray.Dataset - the data to extract the time from
    start_year: int - the start year of the period to extract
    end_year: int - the end year of the period to extract

    Returns:
    dataset: xarray.Dataset - the data with the time extracted
    """
    return dataset.sel(time=slice(str(start_year) + "-12-01", str(end_year) + "-11-30"))


def season_stats(dataset, start_year, end_year, stats=["mean"]):
    """
    Calculate the climatology mean and quantiles for a dataset
    over a given time period for standard meteorological seasons
    and the categories of the terciles for each season.

    Parameters
        dataset: xarray.Dataset - the data to calculate the stats on
        start_year: int - the start year of the period to calculate the stats
        end_year: int - the end year of the period to calculate the stats
        domain: list - the domain to calculate the stats over [lat_min, lat_max, lon_min, lon_max] or string with name of domain
        regrid: string - the regrid method to use, default is 'cons_masked'
        stats: list - the statistics to calculate, default is ['mean']
            options: ['mean', 'seas_ts'], seas_ts is the seasonal time series for all years

    Returns
        ds_stats: dictionary - a dictionary with the calculated stats as xarray.Datasets

    """
    # select the time period
    dataset = sel_season_time(dataset, start_year, end_year)

    # seasonal aggregation - meaning
    # if want other aggreations will need more thought, probably using Iris
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
