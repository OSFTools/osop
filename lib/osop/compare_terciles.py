"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.
"""
# Libraries for working with multi-dimensional arrays
import xarray as xr
import pandas as pd
from osop.util import get_tindex, index
import copy

def mme_products(Services,config,productsfcdir):
    """
    Calls each service percentage file to combine for multi-model ensemble.
    
    Parameters:
    Services (list): List of services to combine
    config (dict): The cofiguraiton parameters for the forecast
    productsfcdir (str): The location for the files to output too and get from
    
    Returns:
    None
    Saves array (x-array) - The multi-model ensemble forecast percentages.
    """
    
    #Remove when happy
    del Services["jma"]
    
    months_1m = list(range(len('{leads_str}'.format(**config))))
    
    #Create a empty list for storage of arrays
    del Services["{origin}".format(**config)]
    datasets_3m = []
    
    for value in months_1m:
        datasets_1m = []
        for origin, systemfc in Services.items():
        #open all services and store
            config_copy = copy.deepcopy(config) #Not sure this is better that just updating it in loop at reupdating it at end outside of loop
            config_copy.update({'origin' : origin, 'systemfc': systemfc})
            if origin == "eccc_can":
                config_copy.update({'origin' : "eccc", 'systemfc': '4'})
            if origin == "eccc_gem5":
                config_copy.update({'origin': "eccc", 'systemfc': '5'})

            local_1m = "{fpath}/{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.imonth_{month}.forecast_percentages.nc".format(
                fpath=productsfcdir,month=value, **config_copy)
            print("this is local_1m", local_1m)
            ds_1m = xr.open_dataset(local_1m)
            datasets_1m.append(ds_1m)
        #Stack and average for percantages
        stacked_1m = xr.concat(datasets_1m, dim='model')
        mme_products_1m = stacked_1m.mean(dim='model')

        #Save out
        mme_fname = "{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.imonth_{month}".format(month=value,
        **config)
        mme_products_1m.to_netcdf(f"{productsfcdir}/{mme_fname}.forecast_percentages.nc")

    for origin, systemfc in Services.items():
        #open all services and store
        config_copy = copy.deepcopy(config) #Not sure this is better that just updating it in loop at reupdating it at end outside of loop
        config_copy.update({'origin' : origin, 'systemfc': systemfc})
        if origin == "eccc_can":
            config_copy.update({'origin' : "eccc", 'systemfc': '4'})
        if origin == "eccc_gem5":
            config_copy.update({'origin': "eccc", 'systemfc': '5'})

        local = "{fpath}/{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.3m.forecast_percentages.nc".format(
        fpath=productsfcdir, **config_copy)

        ds_3m = xr.open_dataset(local)
        datasets_3m.append(ds_3m)

    #Stack and average for percantages
    stacked_3m = xr.concat(datasets_3m, dim='model')
    mme_products_3m = stacked_3m.mean(dim='model')

    #Save out
    mme_fname = "{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}".format(
        **config)
    mme_products_3m.to_netcdf(f"{productsfcdir}/{mme_fname}.3m.forecast_percentages.nc")
    


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

def mask_cat(fcst,terciles):
    """
    Creates a boolean mask for where the forecast value falls.
    Parameters: 
    fcst (x-array): The forecast array
    hcst_terciles (x-array): The hindcast terciles

    Returns: 
    cat (lower,higher,middle) masks (array): The boolean mask
    """

    v = list(terciles.data_vars)[0]
    lo, hi = [terciles[v].sel(category=c) for c in (0, 1)]
    return fcst < lo, fcst > hi, (fcst > lo) & (fcst < hi)



def three_month(forecast_data, hindcast_terciles, products_forecast, forecast_fname ):
    """
    Produces the three month tercile forecast data in the form of 
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
    #Calculate average over the months
    fcst_3m = forecast_data.rolling(forecastMonth=3).mean()
    fcst_3m = fcst_3m.isel(forecastMonth=(fcst_3m['forecastMonth'] == 4))
    #Select for data
    fcst = fcst_3m[list(fcst_3m.data_vars)[0]]
    #Form counts based on catagories 
    lower, higher, middle = mask_cat(fcst, hindcast_terciles)
    total_percentage = xr.Dataset({
        'lower': percentage(lower),
        'higher': percentage(higher),
        'middle': percentage(middle)
    })
    #Save out file for plots 
    total_percentage.to_netcdf(f"{products_forecast}/{forecast_fname}.3m.forecast_percentages.nc")
        
    

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
    # Add valid_time to the xr.Dataset
    start_month = pd.to_datetime(forecast_data.start_date.values).month
    forecast_indexed = forecast_data.assign_coords({"start_month": start_month})

    #Select the data out
    data_var = list(forecast_data.data_vars)[0]

    #Select for each month
    for i in range(forecast_data.sizes['forecastMonth']):
        fc_one_month = forecast_data.isel(forecastMonth=i)[data_var]
        hctt = hindcast_terciles.isel(forecastMonth=i)
        fc_append_name = f'month_{i}'
    #For each month form catagories
        lower, higher, middle = mask_cat(fc_one_month, hctt)
        total_percentage = xr.Dataset({
            'lower': percentage(lower),
            'higher': percentage(higher),
            'middle': percentage(middle)
        })
    #Save out file for plots
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
    # 1 Month
    hcst_terciles_1m = "{fpath}/{origin}_{systemhc}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.1m.tercile_thresholds.nc".format(
        fpath=products_hindcast,**config
    )
    hcst_terciles_1m = xr.open_dataset(hcst_terciles_1m)
    # 3 Month
    hcst_terciles_3m = "{fpath}/{origin}_{systemhc}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.3m.tercile_thresholds.nc".format(
        fpath=products_hindcast,**config
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
    three_month(forecast_data, hcst_terciles_3m,products_forecast, forecast_fname)

    

    


