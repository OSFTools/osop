# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

# Libraries for working with multi-dimensional arrays
"""Functions to compare terciles for forecasts and hindcasts."""

import copy

import pandas as pd
import xarray as xr

from osop.util import get_tindex, index


def update_config(origin, systemfc, config):
    """Create a copy of the config dict to be used for repeated load in of tercile forecasts.

    Parameters
    ----------
    Origin (Str): The service to be loaded
    Systemfc (Str): The service version
    config (Dict): The dictionary to be copied and variated

    Returns
    -------
    Config_copy (Dict): The copy of the dictionary with new updated names
    """
    config_copy = copy.deepcopy(config)
    config_copy.update({"origin": origin, "systemfc": systemfc})
    if origin == "eccc_can":
        config_copy.update({"origin": "eccc", "systemfc": "4"})
    elif origin == "eccc_gem5":
        config_copy.update({"origin": "eccc", "systemfc": "5"})
    return config_copy


def mme_process_forecasts(
    months, suffix, Services, services_values, productsfcdir, config, services_weights
):
    """Load each tercile forecast and combines for MME.

    Parameters
    ----------
    months (int): Set to None or value of month based on leads
    suffix (str): Used for naming between 3months or the imonth lead
    Services (list): List of services to combine
    services_values(dict): list of service values i.e. ecmwf 51
    config (dict): The cofiguraiton parameters for the forecast
    productsfcdir (str): The location for the files to output too and get from
    serives_weights(list): the weight of the service in context of mme

    Returns
    -------
    mme_combined (array): The combined mme forecast array.

    """
    mme_combined = None
    n_members = len(Services)
    for origin, systemfc in services_values.items():
        config_copy = update_config(origin, systemfc, config)
        member_weight = services_weights[origin]

        if suffix == "imonth":
            file_name = "{fpath}/{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.imonth_{month}.forecast_percentages.nc".format(
                fpath=productsfcdir, month=months, **config_copy
            )
        else:
            file_name = "{fpath}/{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.3m.forecast_percentages.nc".format(
                fpath=productsfcdir, **config_copy
            )
        ds = xr.open_dataset(file_name)
        if mme_combined is None:
            mme_combined = xr.zeros_like(ds)
        mme_combined += ds * member_weight
    return mme_combined


def mme_products(Services, config, productsfcdir):
    """Load each tercile forecast and combines for MME.

    Parameters
    ----------
    Services (list): List of services to combine
    config (dict): The cofiguraiton parameters for the forecast
    productsfcdir (str): The location for the files to output too and get from

    Returns
    -------
    None
    Saves array (x-array) - The multi-model ensemble forecast percentages.
    """
    # remove mme from the list that's worked on
    del Services["{origin}".format(**config)]
    # Remove when happy
    del Services["jma"]

    services_values = {
        origin: val[0] if isinstance(val, (list, tuple)) else val
        for origin, val in Services.items()
    }
    services_weights = {
        origin: val[1] if isinstance(val, (list, tuple)) and len(val) > 1 else 1
        for origin, val in Services.items()
    }

    # turn the leads string into array to allow naming of individual months
    months_1m = list(range(len("{leads_str}".format(**config))))
    for month in months_1m:
        mme_products_1m = mme_process_forecasts(
            month,
            "imonth",
            Services,
            services_values,
            productsfcdir,
            config,
            services_weights,
        )
        mme_fname_1m = "{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.imonth_{month}".format(
            month=month, **config
        )
        mme_products_1m.to_netcdf(
            f"{productsfcdir}/{mme_fname_1m}.forecast_percentages.nc"
        )

    mme_products_3m = mme_process_forecasts(
        None, "3m", Services, services_values, productsfcdir, config, services_weights
    )
    mme_fname_3m = "{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}".format(
        **config
    )
    mme_products_3m.to_netcdf(
        f"{productsfcdir}/{mme_fname_3m}.3m.forecast_percentages.nc"
    )


def percentage(array):
    """Take a boolean mask for a forecast dataset and return the percentage of Trues.

    Parameters
    ----------
    array (x-array): The input boolean mask

    Returns
    -------
    array_percentage (x-array): The percentage values
    """
    array_concact = xr.concat(array, dim="number")
    array_sum = array_concact.sum(dim="number")
    array_percentage = (array_sum / array_concact.sizes["number"]) * 100

    return array_percentage


def mask_cat(fcst, terciles):
    """Create a boolean mask for where the forecast value falls.

    Parameters
    ----------
    fcst (x-array): The forecast array
    hcst_terciles (x-array): The hindcast terciles

    Returns
    -------
    cat (lower,higher,middle) masks (array): The boolean mask
    """
    v = list(terciles.data_vars)[0]
    lo, hi = [terciles[v].sel(category=c) for c in (0, 1)]
    return fcst < lo, fcst > hi, (fcst > lo) & (fcst < hi)


def three_month(forecast_data, hindcast_terciles, products_forecast, forecast_fname):
    """Produce a three month tercile forecast.

    This takes data in the form of an xarray that contains the month, the percentage and the lat-lon coordinates.

    Parameters
    ----------
    forecast_data (x-array): The re-indexed forecast data.
    hindcast_terciles (x-array): The x-array that contains the matching tercile categories.
    products_forecast (str): The location for the files to output too.
    forcast_fname (str): The name of the forecast data.

    Returns
    -------
    None
    Saves output data-array that contains the percent values for each tercile and co-ord.

    """
    # Calculate average over the months
    fcst_3m = forecast_data.rolling(forecastMonth=3).mean()
    fcst_3m = fcst_3m.isel(forecastMonth=(fcst_3m["forecastMonth"] == 4))
    # Select for data
    fcst = fcst_3m[list(fcst_3m.data_vars)[0]]

    # Form masks based on categories
    lower, higher, middle = mask_cat(fcst, hindcast_terciles)
    total_percentage = xr.Dataset(
        {
            "lower": percentage(lower),
            "higher": percentage(higher),
            "middle": percentage(middle),
        }
    )
    # Save out file for plots
    total_percentage.to_netcdf(
        f"{products_forecast}/{forecast_fname}.3m.forecast_percentages.nc"
    )


def one_month(forecast_data, hindcast_terciles, products_forecast, forecast_fname):
    """Produce a one month tercile forecast.

    This takes data in the form of an xarray that contains the month, the percentage and the lat-lon coordinates.

    Parameters
    ----------
    forecast_data (x-array): The re-indexed forecast data.
    hindcast_terciles (x-array): The x-array that contains the matching tercile categories.
    products_forecast (str): The location for the files to output too.
    forcast_fname (str): The name of the forecast data.

    Returns
    -------
    None
    Saves output data-array that contains the percent values for each tercile and co-ord.

    """
    # Add valid_time to the xr.Dataset
    start_month = pd.to_datetime(forecast_data.start_date.values).month
    forecast_indexed = forecast_data.assign_coords({"start_month": start_month})

    # Select the data out
    data_var = list(forecast_data.data_vars)[0]

    # Select for each month
    for i in range(forecast_data.sizes["forecastMonth"]):
        fc_one_month = forecast_data.isel(forecastMonth=i)[data_var]
        hctt = hindcast_terciles.isel(forecastMonth=i)
        fc_append_name = f"month_{i}"
        # For each month form categories
        lower, higher, middle = mask_cat(fc_one_month, hctt)
        total_percentage = xr.Dataset(
            {
                "lower": percentage(lower),
                "higher": percentage(higher),
                "middle": percentage(middle),
            }
        )
        # Save out file for plots
        total_percentage.to_netcdf(
            f"{products_forecast}/{forecast_fname}.i{fc_append_name}.forecast_percentages.nc"
        )


def compute_forecast(config, downloaddir, products_hindcast, products_forecast):
    """Calculate tercile forecast data for 1 month forecasts.

    Args:
        config (dict): A dictionary containing the configuration parameters.
        downloaddir (str): The path to the download directory of the forecasts grib.
        products_hindcast (str): The path to the tercile categories from compute_products.
        products_forecast (str): The output location for the products generated.

    Returns
    -------
        None
    """
    # hindcast data set info
    # 1 Month
    hcst_terciles_1m = "{fpath}/{origin}_{systemhc}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.1m.tercile_thresholds.nc".format(
        fpath=products_hindcast, **config
    )
    hcst_terciles_1m = xr.open_dataset(hcst_terciles_1m)
    # 3 Month
    hcst_terciles_3m = "{fpath}/{origin}_{systemhc}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.3m.tercile_thresholds.nc".format(
        fpath=products_hindcast, **config
    )
    hcst_terciles_3m = xr.open_dataset(hcst_terciles_3m)

    # forecast data set info
    forecast_local = "{fpath}/{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.grib".format(
        fpath=downloaddir, **config
    )
    forecast_fname = "{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}".format(
        fpath=downloaddir, **config
    )
    st_dim_name = get_tindex(forecast_local)
    forecast_data = index(forecast_local, st_dim_name)
    one_month(forecast_data, hcst_terciles_1m, products_forecast, forecast_fname)
    three_month(forecast_data, hcst_terciles_3m, products_forecast, forecast_fname)
