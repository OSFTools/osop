# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

# Libraries for working with multi-dimensional arrays
import xarray as xr
import pandas as pd
import numpy as np
import os
import eccodes
import matplotlib.pyplot as plt
from osop.util import get_tindex, index
from osop.compare_terciles import update_config, mme_process_forecasts



# Date and calendar libraries
from dateutil.relativedelta import relativedelta
import logging


def process_mme_products(array, output, aggr, config, sig, member_weight=0.1):
    """
    Calculate anomalies and save them to netCDF files.

    Parameters:
    array (x-array): Target array.
    output (x-array): Enters as an empty array
    aggr (str): 1m or 3m array
    config (dic): Configuration parameters.
    sig (str): Signifier for addition to save path for mean, anom or tercile mme.
    member_weight (float): The fractional weight for the service, i.e. weight/sum of weights
    Returns:
    Returns x-array of the combined mme for 1 month and 3 month combined
    """
    # populate array and divide by members fractional weight
    if output[aggr] is None:
        output[aggr] = xr.zeros_like(array)
    output[aggr] += (array * member_weight)

    save_name = "{origin}_{systemfc}_1993-2016_monthly_mean_{start_month}_{leads_str}_{area_str}_{var}.{aggr}.{sig}".format(
        **config,
        aggr=aggr,
        sig=sig,
    )

    return output[aggr], save_name


def mme_products_hindcast(services, config, productsdir):
    """
    Loads each indexed datasets and combines for mme.

    Parameters:
    Services (list): List of services to combine
    config (dict): The cofiguraiton parameters for the forecast
    productsfcdir (str): The location for the files to output too and get from

    Returns:
    None
    Saves array (x-array) - The multi-model ensemble forecast percentages.
    """
    # remove mme from the list thats worked on
    del services["{origin}".format(**config)]
    # Remove when jma regridded
    del services["jma"]

    mme_combined = {}
    mme_combined_mean = {}
    mme_combined_anom = {}


    services_values = {origin: val[0] if isinstance(val, (list, tuple)) else val
                        for origin, val in services.items()}
    services_weights = {origin: val[1] if isinstance(val, (list, tuple)) and len(val) > 1 else 1
                        for origin, val in services.items()}

    
    for aggr in ["1m", "3m"]:
        mme_combined[aggr] = None
        mme_combined_mean[aggr] = None
        mme_combined_anom[aggr] = None

        targets = {
            "tercile_probs.nc": mme_combined,
            "ensmean.nc": mme_combined_mean,
            "anom.nc": mme_combined_anom,
        }
        for suffix, target in targets.items():
            for origin, system_value in services_values.items():
                config_copy_hc = update_config(origin, system_value, config)

                # Load and run each array, 1m-3m and tercile, mean and anom
                file_name = "{fpath}/{origin}_{systemfc}_1993-2016_monthly_mean_{start_month}_{leads_str}_{area_str}_{var}.{aggr}.{suffix}".format(
                    fpath=productsdir,
                    **config_copy_hc,
                    aggr=aggr,
                    suffix=suffix,
                )
                logging.debug(
                    "Combining MME product for {suffix}".format(suffix=suffix)
                )
                output, save_name = process_mme_products(
                    xr.open_dataset(file_name),
                    target,
                    aggr,
                    config,
                    suffix,
                    member_weight=services_weights[origin]
                )
            # Save out final mme array
            output.to_netcdf(f"{productsdir}/{save_name}")
    return


def calc_anoms(hcst, hcst_bname, config, productsdir):
    """
    Calculate anomalies and save them to netCDF files.

    Parameters:
    hcst_fname (str): File path of the hindcast grib file.
    hcst_bname (str): Base name of the hindcast grib file.
    config (dict): Configuration parameters.
    st_dim_name (str): Name of the start date dimension (important for lagged models)
    productsdir (str): Directory path to save the netCDF files.

    Returns:
    saves 1 month and 3 month anomalies to netCDF files
    returns a tuple of xarray datasets containing the original hindcast data and the 3-month aggregated data.
    """
    logging.debug("Re-arranging time metadata in xr.Dataset object")
    # Add start_month to the xr.Dataset
    start_month = pd.to_datetime(hcst.start_date.values[0]).month
    hcst = hcst.assign_coords({"start_month": start_month})

    # Add valid_time to the xr.Dataset
    vt = xr.DataArray(
        dims=("start_date", "forecastMonth"),
        coords={"forecastMonth": hcst.forecastMonth, "start_date": hcst.start_date},
    )
    vt.data = [
        [
            pd.to_datetime(std) + relativedelta(months=fcmonth - 1)
            for fcmonth in vt.forecastMonth.values
        ]
        for std in vt.start_date.values
    ]
    hcst = hcst.assign_coords(valid_time=vt)

    # CALCULATE 3-month AGGREGATIONS
    # NOTE rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
    logging.debug("Computing 3-month aggregation")
    # rollng method defaults to look backwards
    hcst_3m = hcst.rolling(forecastMonth=3).mean()
    # Want only 3 month mean with complete 3 months
    hcst_3m = hcst_3m.where(hcst_3m.forecastMonth >= int(config["leads"][2]), drop=True)

    # Calculate Anomalies (and save to file)
    logging.debug("Computing anomalies 1m")
    hcmean = hcst.mean(["number", "start_date"])
    # Calculate Mean across all ensemble members
    hc_ens_mean = hcst.mean(["number"])
    anom = hcst - hcmean
    anom = anom.assign_attrs(reference_period="{hcstarty}-{hcendy}".format(**config))

    logging.debug("Computing anomalies 3m")
    hcmean_3m = hcst_3m.mean(["number", "start_date"])
    hc_ens_mean_3m = hcst_3m.mean(["number"])
    anom_3m = hcst_3m - hcmean_3m
    anom_3m = anom_3m.assign_attrs(
        reference_period="{hcstarty}-{hcendy}".format(**config)
    )

    logging.debug("Saving mean and anomalies 1m/3m to netCDF files")
    anom.to_netcdf(f"{productsdir}/{hcst_bname}.1m.anom.nc")
    anom_3m.to_netcdf(f"{productsdir}/{hcst_bname}.3m.anom.nc")
    hc_ens_mean.to_netcdf(f"{productsdir}/{hcst_bname}.1m.ensmean.nc")
    hc_ens_mean_3m.to_netcdf(f"{productsdir}/{hcst_bname}.3m.ensmean.nc")

    return hcst, hcst_3m


def get_thresh(icat, quantiles, xrds, dims=["number", "start_date"]):
    """Calculate the boundaries of forecast categories defined by quantiles e.g. terciles

    Args:
        icat (int): The category number. 0 (lower than), 1 (normal), 2 (higher than).
        quantiles (list): The list of quantiles. Use [1/3., 2/3.] for terciles.
        xrds (xarray.Dataset): The dataset containing the hindcast data.
        dims (list, optional): The dimensions to consider when calculating the quantiles.
            Defaults to ['number', 'start_date'].

    Returns:
        tuple: A tuple containing the lower and upper boundaries for the forecast category.
    """

    if not all(elem in xrds.sizes for elem in dims):
        raise ValueError(
            "Some of the dimensions in {} is not present in the xr.Dataset {}".format(
                dims, xrds
            )
        )
    else:
        if icat == 0:
            xrds_lo = -np.inf
            xrds_hi = xrds.quantile(quantiles[icat], dim=dims)

        elif icat == len(quantiles):
            xrds_lo = xrds.quantile(quantiles[icat - 1], dim=dims)
            xrds_hi = np.inf

        else:
            xrds_lo = xrds.quantile(quantiles[icat - 1], dim=dims)
            xrds_hi = xrds.quantile(quantiles[icat], dim=dims)

    return xrds_lo, xrds_hi


def prob_terc(config, hcst_bname, hcst, hcst_3m, productsdir):
    """
    Calculate probabilities for tercile categories
    by counting members within each category and save them to netCDF files.
    This function computes the tercile thresholds for both 1-month and 3-month aggregated hindcast data.
    and saves them to netCDF files.

    Args:
        hcst_bname (str): Basename of hincast file.
        hcst (xarray.Dataset): The dataset containing the hindcast data.
        hcst_3m (xarray.Dataset): The dataset containing the 3-month aggregated hindcast data.
        productsdir (str): Directory path to save the netCDF files.

    Returns:
        NA
        Saves tercile forecasts to netcdf file.
    """

    logging.debug("Computing probabilities (tercile categories)")
    quantiles = [1 / 3.0, 2 / 3.0]
    numcategories = len(quantiles) + 1

    for aggr, h in [("1m", hcst), ("3m", hcst_3m)]:
        if os.path.isfile(f"{productsdir}/{hcst_bname}.{aggr}.tercile_probs.nc"):
            logging.debug(f"{productsdir}/{hcst_bname}.{aggr}.tercile_probs.nc exists")
        else:
            logging.debug(f"Computing tercile probabilities {aggr}")

            l_probs_hcst = list()
            # Loop over the categories and calculate the probabilities
            for icat in range(numcategories):

                h_lo, h_hi = get_thresh(icat, quantiles, h)
                probh = np.logical_and(h > h_lo, h <= h_hi).sum("number") / float(
                    h.sizes["number"]
                )

                # Instead of using the coordinate 'quantile' coming from the hindcast xr.Dataset
                # we will create a new coordinate called 'category'
                if "quantile" in probh:
                    probh = probh.drop("quantile")
                l_probs_hcst.append(probh.assign_coords({"category": icat}))

                # on second iteration the values of h_lo and h_hi are the
                # quantiles we wish to save
                if icat == 1:
                    tercs = xr.concat([h_lo, h_hi], dim="category")
                    if "quantile" in tercs:
                        tercs = tercs.drop("quantile")
                    # include metadata about the reference period and start month
                    tercs = tercs.assign_attrs(
                        reference_period="{hcstarty}-{hcendy}".format(**config)
                    )
                    tercs = tercs.assign_attrs(start_month=f"{config['start_month']}")

            logging.debug(f"Concatenating {aggr} tercile probs categories")
            probs = xr.concat(l_probs_hcst, dim="category")
            logging.debug(f"Saving {aggr} tercile probs netCDF files")
            probs.to_netcdf(f"{productsdir}/{hcst_bname}.{aggr}.tercile_probs.nc")

            logging.debug(f"Saving tercile thresholds {aggr} netCDF files")
            tercs.to_netcdf(f"{productsdir}/{hcst_bname}.{aggr}.tercile_thresholds.nc")


def calc_products(config, downloaddir, productsdir):
    """Calculate anomalies and tercile probabilities for a given hindcast dataset

    Args:
        config (dict): Configuration parameters.
        downloaddir (str): Directory path to save the netCDF files.
    """

    hcst_bname = "{origin}_{system}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{var}".format(
        **config
    )
    hcst_fname = f"{downloaddir}/{hcst_bname}.grib"

    # For the re-shaping of time coordinates in xarray.Dataset we need to select the right coord
    #  -> burst mode ensembles (e.g. ECMWF SEAS5) use "time". This is the default option
    #  -> lagged start ensembles (e.g. MetOffice GloSea6) use "indexing_time" (see CDS documentation about nominal start date)
    st_dim_name = get_tindex(hcst_fname)
    hcst = index(hcst_fname, st_dim_name)
    logging.debug("this is hcst from index", hcst)

    ## calc anoms
    logging.info(f"Calculating anomalies for {hcst_bname}")
    hcst, hcst_3m = calc_anoms(hcst, hcst_bname, config, productsdir)
    ## calc terc probs and thresholds
    logging.info(f"Calculating tercile probabilities for {hcst_bname}")
    prob_terc(config, hcst_bname, hcst, hcst_3m, productsdir)


def calc_products_mme(services, config, productsdir):
    """
    Calculate anomalies and tercile probabilities for a the mme combined data.

    Args:
        Services (dict): All the services intended for the mme.
        config (dict): Configuration parameters.
        downloaddir (str): Directory path to save the netCDF files.
    """

    hcst_bname = "{origin}_{systemfc}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{var}".format(
        **config
    )

    # Combined the services to produce a mme for hindcast verifcation
    hcst = mme_products_hindcast(services, config, productsdir)
