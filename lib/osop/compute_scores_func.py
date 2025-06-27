"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.
"""

# This script currently can only be used to create scores using ERA5 as the comparison dataset

# Libraries for working with multi-dimensional arrays
import xarray as xr
import pandas as pd
import numpy as np


# Forecast verification metrics with xarray
import xskillscore as xs
# Regridding packages needed for JMA
import xesmf as xe

# import needed local functions
from osop.compute_products_func import get_thresh


def read_obs(obs_fname, config):
    """
    Read observation data from a file and calculate 3 month averages.

    Parameters:
    obs_fname (str): The file path of the observation data file.
    config (dict): A dictionary containing configuration parameters.

    Returns:
    obs_ds (xarray.Dataset): Preprocessed observation data with monthly time resolution.
    obs_ds_3m (xarray.Dataset): Preprocessed observation data with 3-monthly time resolution.
    """

    obs_ds = xr.open_dataset(obs_fname, engine="cfgrib")

    # Renaming to match hindcast names
    obs_ds = obs_ds.rename(
        {"latitude": "lat", "longitude": "lon", "time": "start_date"}
    ).swap_dims({"start_date": "valid_time"})

    # Assign 'forecastMonth' coordinate values
    fcmonths = [
        mm + 1 if mm >= 0 else mm + 13
        for mm in [
            t.month - config["start_month"]
            for t in pd.to_datetime(obs_ds.valid_time.values)
        ]
    ]
    obs_ds = obs_ds.assign_coords(forecastMonth=("valid_time", fcmonths))
    # Drop obs values not needed (earlier than first start date)
    # to create well shaped 3-month aggregations from obs.
    obs_ds = obs_ds.where(
        obs_ds.valid_time
        >= np.datetime64("{hcstarty}-{start_month:02d}-01".format(**config)),
        drop=True,
    )

    # Calculate 3-month aggregations
    # NOTE rolling() assigns the label to the end of the N month period
    print("Calculate observation 3-monthly aggregations")
    obs_ds_3m = obs_ds.rolling(valid_time=3).mean()
    obs_ds_3m = obs_ds_3m.where(
        obs_ds_3m.forecastMonth >= int(config["leads"][2]), drop=True
    )

    # As we don't need it anymore, we can safely remove 'forecastMonth'
    obs_ds = obs_ds.drop("forecastMonth")
    obs_ds_3m = obs_ds_3m.drop("forecastMonth")

    obs_ds["valid_time"] = obs_ds.valid_time.dt.strftime("%Y-%m")
    obs_ds_3m["valid_time"] = obs_ds_3m.valid_time.dt.strftime("%Y-%m")
    return obs_ds, obs_ds_3m


def swap_dims(hcst, obs):
    """
    Swaps dimensions appropriately and alligns hindcast with observed dataset
    for deterministic verification measures.

    Parameters:
    hcst: the hindcast dataset.
    obs: the matching observed dataset.

    Returns:
    thishcst: the alligned hindcast dataset.
    thisobs: the alligned observed dataset.
    """

    for this_fcmonth in hcst.forecastMonth.values:
        print(f"forecastMonth={this_fcmonth}")
        thishcst = hcst.sel(forecastMonth=this_fcmonth).swap_dims(
            {"start_date": "valid_time"}
        )
        thishcst["valid_time"] = thishcst.valid_time.dt.strftime("%Y-%m")
        print(thishcst["valid_time"])
        print(obs["valid_time"])
        thisobs = obs.where(obs.valid_time == thishcst.valid_time, drop=True)

        return thishcst, thisobs

def align_data(thishcst, this_obs):
    """
    Regrids the dataset appropriatley for its type (planned expansion for percipertation)

    Parameters:
    thishcst (xarray.Dataset): Hindcast data for regrid
    thisobs (xarray.Dataset): Matching obs dataset to be set as the target

    Returns:
    thishcst (xarray.Dataset): Regridded dataset to be used for analysis
    thisobs (xarray.Dataset): Matching obs dataset (no changes)
    """
    try:
        regridder =  xe.Regridder(thishcst, this_obs, "bilinear")
        thishcst = regridder(thishcst, keep_attrs=True)
    except Exception as e:
                print(f"Alignment failed {e}: {e}")
                raise KeyError("Alignment failed: please check dataset entry")    
    
    return thishcst, this_obs



def scores_dtrmnstc(obs_ds, obs_ds_3m, hcst_bname, scoresdir, productsdir):
    """
    Compute deterministic scores.

    Parameters:
    obs_ds (xarray.Dataset): Observation / reanalysis data, monthly resolution.
    obs_ds_3m (xarray.Dataset): Observation / reanalysis 3-month aggregated data.
    hcst_bname (str): Basename of the hindcast data.
    scoresdir (str): Directory to save the output files.
    productsdir (str): Directory to fetch files from.

    Returns:
    None
    Saves spearman and pearson correlation and p-value to netCDF files.
    """

    # Loop over aggregations
    for aggr in ["1m", "3m"]:

        if aggr == "1m":
            o = obs_ds
        elif aggr == "3m":
            o = obs_ds_3m
        else:
            raise ValueError(f"Unknown aggregation {aggr}")

        print(f"Computing deterministic scores for {aggr}-aggregation")

        # Reading mean file
        h = xr.open_dataset(f"{productsdir}/{hcst_bname}.{aggr}.ensmean.nc")
        # Reading anomalies file
        ha = xr.open_dataset(f"{productsdir}/{hcst_bname}.{aggr}.anom.nc")
        is_fullensemble = "number" in ha.dims

        # create empty list to store correlations and p-values to be concatenated after looping over months
        l_corr = list()
        l_corr_pval = list()
        r_corr = list()
        r_corr_pval = list()

        # match dimensions and calculate obs anomalies
        thishcst_em_mean, this_obs_m_match = swap_dims(h, o)
        thishcst_em_anom, this_obs_a_match = swap_dims(ha, o)
        obsmean = this_obs_a_match.mean()
        this_obs_anom = this_obs_a_match - obsmean
        thishcst_em_anom = (
            thishcst_em_anom if not is_fullensemble else thishcst_em_anom.mean("number")
        )
         # Regrid if lattitude or longitude on a varied resolution or grid. 
        if not thishcst_em_mean['lat'].equals(this_obs_m_match['lat']) or not thishcst_em_mean['lon'].equals(this_obs_m_match['lon']):
            thishcst_em_mean, this_obs_m_match = align_data(thishcst_em_mean, this_obs_m_match)
        if not thishcst_em_anom['lat'].equals(this_obs_anom['lat']) or not thishcst_em_anom['lon'].equals(this_obs_anom['lon']):
            thishcst_em_anom, this_obs_anom = align_data(thishcst_em_anom, this_obs_anom)


        # calculate measures
        l_corr.append(
            xs.spearman_r(thishcst_em_mean, this_obs_m_match, dim="valid_time")
        )
        l_corr_pval.append(
            xs.spearman_r_p_value(thishcst_em_mean, this_obs_m_match, dim="valid_time")
        )
        r_corr.append(xs.pearson_r(thishcst_em_anom, this_obs_anom, dim="valid_time"))
        r_corr_pval.append(
            xs.pearson_r_p_value(thishcst_em_anom, this_obs_anom, dim="valid_time")
        )

        # concatenate correlations and p-values
        corr = xr.concat(l_corr, dim="forecastMonth")
        corr_pval = xr.concat(l_corr_pval, dim="forecastMonth")
        r_corr = xr.concat(r_corr, dim="forecastMonth")
        r_corr_pval = xr.concat(r_corr_pval, dim="forecastMonth")

        print(f"Saving to netCDF file correlation for {aggr}-aggregation")
        corr.to_netcdf(f"{scoresdir}/{hcst_bname}.{aggr}.spearman_corr.nc")
        corr_pval.to_netcdf(
            f"{scoresdir}/{hcst_bname}.{aggr}.spearman_corr_pval.nc"
        )
        r_corr.to_netcdf(f"{scoresdir}/{hcst_bname}.{aggr}.pearson_corr.nc")
        r_corr_pval.to_netcdf(
            f"{scoresdir}/{hcst_bname}.{aggr}.pearson_corr_pval.nc"
        )


def scores_prblstc(obs_ds, obs_ds_3m, hcst_bname, scoresdir, productsdir):
    """
    Compute probabilistic scores and save the results to NetCDF files.

    Parameters:
    obs_ds(xarray.Dataset): Observation / Reanalysis monthly data.
    obs_ds_3m (xarray.Dataset): Observation / Reanalysis 3-month aggregated data.
    hcst_bname (str): Basename of the hindcast probabilities file.
    scoresdir (str): Directory to save the output NetCDF files.
    productsdir (str): Directory to fetch input files from.

    Returns:
    None
    Saves probabilistic scores to NetCDF files.
    """
    # Define quantiles for tercile categories
    quantiles = [1 / 3.0, 2 / 3.0]
    numcategories = len(quantiles) + 1
    # Loop over aggregations
    for aggr in ["1m", "3m"]:
        if aggr == "1m":
            o = obs_ds
        elif aggr == "3m":
            o = obs_ds_3m
        else:
            raise BaseException(f"Unknown aggregation {aggr}")

        print(f"Computing probabilistic scores for {aggr}-aggregation")

        # Read hindcast probabilities file
        probs_hcst = xr.open_dataset(
            f"{productsdir}/{hcst_bname}.{aggr}.tercile_probs.nc"
        )

        l_roc = list()
        l_rps = list()
        l_rocss = list()
        l_bs = list()
        l_rel = list()

        for this_fcmonth in probs_hcst.forecastMonth.values:
            thishcst = probs_hcst.sel(forecastMonth=this_fcmonth).swap_dims(
                {"start_date": "valid_time"}
            )

            # Calculate tercile categories from observations
            l_probs_obs = list()

            # format valid_time to match between obs and hcst
            thishcst["valid_time"] = thishcst.valid_time.dt.strftime("%Y-%m")

            # Select relevant observations
            thiso = o.where(o.valid_time == thishcst.valid_time, drop=True)

            for icat in range(numcategories):

                o_lo, o_hi = get_thresh(icat, quantiles, thiso, dims=["valid_time"])
                probo = 1.0 * np.logical_and(thiso > o_lo, thiso <= o_hi)
                if "quantile" in probo:
                    probo = probo.drop("quantile")
                l_probs_obs.append(probo.assign_coords({"category": icat}))

            thisobs = xr.concat(l_probs_obs, dim="category")

            # Regrid if lattitude or longitude on a varied resolution or grid. 
            if not thishcst['lat'].equals(thisobs['lat']) or not thishcst['lon'].equals(thisobs['lon']):
               thishcst, thisobs = align_data(thishcst,thisobs)

            # Calculate the probabilistic scores
            thisroc = xr.Dataset()
            thisrps = xr.Dataset()
            thisrocss = xr.Dataset()
            thisbs = xr.Dataset()
            thisrel = xr.Dataset()
            for var in thishcst.data_vars:
                var_obs = var

                thisroc[var] = xs.roc(
                    thisobs[var_obs],
                    thishcst[var],
                    dim="valid_time",
                    bin_edges=np.linspace(0, 1, 101),
                )
                thisrocss[var] = (thisroc[var] - 0.5) / (1.0 - 0.5)

                thisrps[var] = xs.rps(
                    thisobs[var_obs],
                    thishcst[var],
                    dim="valid_time",
                    category_edges=None,
                    input_distributions="p",
                )

                bscat = list()
                relcat = list()
                for cat in thisobs[var_obs].category:

                    thisobscat = thisobs[var_obs].sel(category=cat)
                    thishcstcat = thishcst[var].sel(category=cat)

                    bscat.append(
                        xs.brier_score(thisobscat, thishcstcat, dim="valid_time")
                    )
                    relcat.append(xs.reliability(thisobscat.astype(int), thishcstcat))

                thisbs[var] = xr.concat(bscat, dim="category")
                thisrel[var] = xr.concat(relcat, dim="category")

            l_roc.append(thisroc)
            l_rps.append(thisrps)
            l_rocss.append(thisrocss)
            l_bs.append(thisbs)
            l_rel.append(thisrel)
            ## Only want first month if aggr == 1m
            if aggr == "1m":
                break

        roc = xr.concat(l_roc, dim="forecastMonth")
        rps = xr.concat(l_rps, dim="forecastMonth")
        rocss = xr.concat(l_rocss, dim="forecastMonth")
        bs = xr.concat(l_bs, dim="forecastMonth")
        rel = xr.concat(l_rel, dim="forecastMonth")

        print("writing scores to netcdf")
        rps.to_netcdf(f"{scoresdir}/{hcst_bname}.{aggr}.rps.nc")
        bs.to_netcdf(f"{scoresdir}/{hcst_bname}.{aggr}.bs.nc")
        roc.to_netcdf(f"{scoresdir}/{hcst_bname}.{aggr}.roc.nc")
        rocss.to_netcdf(f"{scoresdir}/{hcst_bname}.{aggr}.rocss.nc")
        rel.to_netcdf(f"{scoresdir}/{hcst_bname}.{aggr}.rel.nc")


def calc_scores(config, downloaddir, scoresdir, productsdir):
    """
    Calls code to calculate deterministic and probabilistic verification scores.

    Args:
        config (dict): A dictionary containing the configuration parameters.
        downloaddir (str): The path to the download directory.

    Returns:
        None
    """

    hcst_bname = "{origin}_{system}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}".format(
        **config
    )

    ## obs info
    obs_fname = "{fpath}/{obs_name}_{hc_var}_{hcstarty}-{hcendy}_monthly_{start_month}_{obs_str}_{area_str}.grib".format(
        fpath=downloaddir, **config
    )

    ## read obs
    obs_ds, obs_ds_3m = read_obs(obs_fname, config)

    # if the observations are ERA5, rename tprate
    # and convert from m/day to m/s
    if config["obs_name"] == "era5" and config["hc_var"] == "total_precipitation":
        obs_ds = obs_ds.rename({"tp": "tprate"})
        obs_ds_3m = obs_ds_3m.rename({"tp": "tprate"})
        obs_ds["tprate"].attrs["units"] = "m/s"
        obs_ds["tprate"] = obs_ds["tprate"] * 3600 * 24
        obs_ds_3m["tprate"] = obs_ds_3m["tprate"] * 3600 * 24
        obs_ds_3m["tprate"].attrs["units"] = "m/s"
    ## calc scores
    scores_dtrmnstc(obs_ds, obs_ds_3m, hcst_bname, scoresdir, productsdir)
    scores_prblstc(obs_ds, obs_ds_3m, hcst_bname, scoresdir, productsdir)
