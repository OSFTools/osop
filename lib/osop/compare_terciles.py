"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.
"""
# Libraries for working with multi-dimensional arrays
import xarray as xr
import xarray as xr
import pandas as pd
import eccodes

def get_tindex(infile):
    """
    Use eccodes to check if there is an indexing time dimension

    Input
        infile(str): name of file to check
    Returns:
        st_dim_name (str): name of time dimension to use for indexing.
            time for burst ensmeble and indexing_time for lagged

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

def index(forecast_local, st_dim_name):
    """
    Reindex and restyle the forcast grib so that the data layout is consistent
    and compatiable with hindcast terciles. 

    Parameters:
    forecast_local (str): File location for the grib file.
    st_dim_name (str): Name of the start date dimension (important for lagged models)

    Returns:
    A re-indexed x-array for forecast data. 
    """

    print("Reading Forecast data from file")
    forecast_data = xr.open_dataset(
        forecast_local,
        engine="cfgrib",
        backend_kwargs=dict(time_dims=("forecastMonth", st_dim_name)),
    )
    # force dask.array using chunks on leadtime, latitude and longitude coordinate
    forecast_data = forecast_data.chunk({"forecastMonth": 1, "latitude": "auto", "longitude": "auto"})
    forecast_data = forecast_data.rename(
        {"latitude": "lat", "longitude": "lon", st_dim_name: "start_date"}
    )

    print("Re-arranging time metadata in xr.Dataset object")
    # Add start_month to the xr.Dataset
    start_month = pd.to_datetime(forecast_data.start_date.values).month
    forecast_indexed = forecast_data.assign_coords({"start_month": start_month})
    # Add valid_time to the xr.Dataset
    return(forecast_indexed)

def percentage(array):
    """
    Takes a boolean mask for a forecast dataset and returns the percentage of Trues

    Parameters: 
    array (x-array): The input boolean mask

    Returns: 
    array_percentage (x-array): The percentage values
    """
    array_concact = xr.concat(array, dim='number')
    array_sum = array_concact.sum(dim='number')
    array_percentage = (array_sum / array_concact.sizes['number']) * 100

    return array_percentage
    

def one_month(forecast_data, hindcast_terciles, products_forecast, forecast_fname):
    """
    Produces the one month tercile forecast data in the form of 
    an x-array that contains the month, the percentage and the lat-lon coordinates. 

    Parameters:
    forecast_data (x-array): The re-indexed forecast data.
    hindcast_terciles (x-array): The x-array that contains the matching tercile catagories.
    products_forecast (str): The location for the files to output too. 
    forcast_fname (str): The name of the forecast data. 

    Returns:
    None 
    Saves output data-array that contains the percent values for each tercile and co-ord. 
    
    """
    for i in range(forecast_data.sizes['forecastMonth']):
        #Select for each forecast Month
        fc_one_month = forecast_data.isel(forecastMonth=i)
        hctt = hindcast_terciles.isel(forecastMonth=i)

        fc_append_name = f'month_{i}'

        
        #Select for variable
        data_var = list(fc_one_month.data_vars)[0]
        fc_one_month = fc_one_month[data_var]

        data_var = list(hctt.data_vars)[0]

        #Define Catagories
        cat_lower_thresh = hctt[data_var].sel(category=0)
        cat_higher_thresh = hctt[data_var].sel(category=1)

        #Create empty lists to store into
        cat_lower_masks = []
        cat_higher_masks = []
        cat_middle_masks = []

        #for each member - assign catagory
        for i in range(fc_one_month.sizes['number']):
            fc_one_month_slice = fc_one_month.isel(number=i)
            mask = fc_one_month_slice < cat_lower_thresh
            cat_lower_masks.append(mask)
            mask1 = fc_one_month_slice > cat_higher_thresh
            cat_higher_masks.append(mask1)
            mask3 = (fc_one_month_slice > cat_lower_thresh) & (fc_one_month_slice < cat_higher_thresh)
            cat_middle_masks.append(mask3)
    
        #"Percentagise" the boolean masks
        percentage_lower= percentage(cat_lower_masks)
        percentage_higher = percentage(cat_higher_masks)
        percentage_middle = percentage(cat_middle_masks)

        #Reorganise to one dataset for ease
        total_percentage = xr.Dataset({
            'lower':percentage_lower,
            'higher':percentage_higher,
            'middle':percentage_middle
        })

        total_percentage.to_netcdf(f"{products_forecast}/{forecast_fname}.i{fc_append_name}.forecast_percentages.nc")


def compute_forecast(config, downloaddir, products_hindcast, products_forecast):
    """
    Calls code to calculate tercile forecast data for 1 month forecasts.

    Args:
        config (dict): A dictionary containing the configuration parameters.
        downloaddir (str): The path to the download directory of the forecasts grib.
        products_hindcast (str): The path to the tercile catagories from compute_products.
        products_forecast (str): The output location for the products generated.

    Returns:
        None
    """
    # hindcast data set info 
    hcst_terciles = "{fpath}/{origin}_{systemhc}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.1m.tercile_thresholds.nc".format(
        fpath=products_hindcast,**config
    )
    hcst_terciles = xr.open_dataset(hcst_terciles)

    # forecast data set info
    forecast_local = "{fpath}/{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.grib".format(
        fpath=downloaddir, **config
    )
    forecast_fname = "{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}".format(
        fpath=downloaddir, **config
    )
    st_dim_name = get_tindex(forecast_local)
    forecast_data = index(forecast_local, st_dim_name)
    one_month(forecast_data, hcst_terciles, products_forecast, forecast_fname)

    

    


