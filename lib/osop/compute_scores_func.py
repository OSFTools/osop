# This script currently can only be used to create scores using ERA5 as the comparison dataset

# Libraries for working with multi-dimensional arrays
import xarray as xr
import pandas as pd
import numpy as np


# Forecast verification metrics with xarray
import xskillscore as xs

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


def scores_dtrmnstc(obs_ds, obs_ds_3m, hcst_bname, downloaddir):
    """
    Compute deterministic scores.

    Parameters:
    obs_ds (xarray.Dataset): Observation / reanalysis data, monthly resolution.
    obs_ds_3m (xarray.Dataset): Observation / reanalysis 3-month aggregated data.
    hcst_bname (str): Basename of the hindcast data.
    downloaddir (str): Directory to save the output files.

    Returns:
    None
    Saves spearman correlation and p-value to netCDF files.
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

        # Read anomalies file
        h = xr.open_dataset(f"{downloaddir}/{hcst_bname}.{aggr}.anom.nc")
        is_fullensemble = "number" in h.dims

        # create empty list to store correlations and p-values to be concatenated after looping over months
        l_corr = list()
        l_corr_pval = list()

        for this_fcmonth in h.forecastMonth.values:
            print(f"forecastMonth={this_fcmonth}")
            thishcst = h.sel(forecastMonth=this_fcmonth).swap_dims(
                {"start_date": "valid_time"}
            )
            thishcst["valid_time"] = thishcst.valid_time.dt.strftime("%Y-%m")
            print(thishcst["valid_time"])
            print(o["valid_time"])
            thisobs = o.where(o.valid_time == thishcst.valid_time, drop=True)
            thishcst_em = thishcst if not is_fullensemble else thishcst.mean("number")
            l_corr.append(xs.spearman_r(thishcst_em, thisobs, dim="valid_time"))
            l_corr_pval.append(
                xs.spearman_r_p_value(thishcst_em, thisobs, dim="valid_time")
            )

        # concatenate correlations and p-values
        corr = xr.concat(l_corr, dim="forecastMonth")
        corr_pval = xr.concat(l_corr_pval, dim="forecastMonth")

        print(f"Saving to netCDF file correlation for {aggr}-aggregation")
        corr.to_netcdf(f"{downloaddir}/scores/{hcst_bname}.{aggr}.spearman_corr.nc")
        corr_pval.to_netcdf(
            f"{downloaddir}/scores/{hcst_bname}.{aggr}.spearman_corr_pval.nc"
        )


def scores_prblstc(obs_ds, obs_ds_3m, hcst_bname, downloaddir):
    """
    Compute probabilistic scores and save the results to NetCDF files.

    Parameters:
    obs_ds(xarray.Dataset): Observation / Reanalysis monthly data.
    obs_ds_3m (xarray.Dataset): Observation / Reanalysis 3-month aggregated data.
    hcst_bname (str): Basename of the hindcast probabilities file.
    downloaddir (str): Directory to save the output NetCDF files.

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
            f"{downloaddir}/{hcst_bname}.{aggr}.tercile_probs.nc"
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
        rps.to_netcdf(f"{downloaddir}/scores/{hcst_bname}.{aggr}.rps.nc")
        bs.to_netcdf(f"{downloaddir}/scores/{hcst_bname}.{aggr}.bs.nc")
        roc.to_netcdf(f"{downloaddir}/scores/{hcst_bname}.{aggr}.roc.nc")
        rocss.to_netcdf(f"{downloaddir}/scores/{hcst_bname}.{aggr}.rocss.nc")
        rel.to_netcdf(f"{downloaddir}/scores/{hcst_bname}.{aggr}.rel.nc")


def calc_scores(config, downloaddir):
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
    scores_dtrmnstc(obs_ds, obs_ds_3m, hcst_bname, downloaddir)
    scores_prblstc(obs_ds, obs_ds_3m, hcst_bname, downloaddir)
