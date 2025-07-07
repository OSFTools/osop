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

def index(hcst_fname, st_dim_name):
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

    ##print("Reading HCST data from file")
    hcst = xr.open_dataset(
        hcst_fname,
        engine="cfgrib",
        backend_kwargs=dict(time_dims=("forecastMonth", st_dim_name)),
    )
    # force dask.array using chunks on leadtime, latitude and longitude coordinate
    hcst = hcst.chunk({"forecastMonth": 1, "latitude": "auto", "longitude": "auto"})
    hcst = hcst.rename(
        {"latitude": "lat", "longitude": "lon", st_dim_name: "start_date"}
    )

    ##print("Re-arranging time metadata in xr.Dataset object")
    # Add start_month to the xr.Dataset
    start_month = pd.to_datetime(hcst.start_date.values).month
    hcst = hcst.assign_coords({"start_month": start_month})
    # Add valid_time to the xr.Dataset
    return(hcst)


def one_month(forecast_data, hindcast_terciles, products_forecast, forecast_fname):
    for i in range(forecast_data.sizes['forecastMonth']):
        #Select for each forecast Month
        fc = forecast_data.isel(forecastMonth=i)
        hctt = hindcast_terciles.isel(forecastMonth=i)

        

        #Select for variable
        data_var = list(fc.data_vars)[0]
        t2mfc = fc[data_var]

        data_var = list(hctt.data_vars)[0]

        #Define Catagories
        cat_lower_thresh = hctt[data_var].sel(category=0)
        cat_higher_thresh = hctt[data_var].sel(category=1)

        #Create empty lists to store into
        cat_lower_masks = []
        cat_higher_masks = []
        cat_middle_masks = []

        #for each member - assign catagory
        for i in range(t2mfc.sizes['number']):
            t2m_slice = t2mfc.isel(number=i)
            mask = t2m_slice < cat_lower_thresh
            cat_lower_masks.append(mask)
            mask1 = t2m_slice > cat_higher_thresh
            cat_higher_masks.append(mask1)
            mask3 = (t2m_slice > cat_lower_thresh) & (t2m_slice < cat_higher_thresh)
            cat_middle_masks.append(mask3)
    
        #reassmeble
        cat_lower_mask_array = xr.concat(cat_lower_masks, dim='number')
        cat_higher_mask_array = xr.concat(cat_higher_masks, dim='number')
        cat_middle_mask_array = xr.concat(cat_middle_masks, dim='number')

        #create counts
        true_counts_lower = cat_lower_mask_array.sum(dim='number')
        true_counts_higher = cat_higher_mask_array.sum(dim='number')
        true_counts_middle = cat_middle_mask_array.sum(dim='number')

        #create percentages
        percentage_true_lower = (true_counts_lower / cat_lower_mask_array.sizes['number']) * 100
        percentage_true_higher = (true_counts_higher / cat_higher_mask_array.sizes['number']) * 100
        percentage_true_middle = (true_counts_middle / cat_middle_mask_array.sizes['number']) * 100

        print("this is percentage_true_lower:",percentage_true_lower)
        percentage_true_lower.to_netcdf(f"{products_forecast}/{forecast_fname}.forecast_lower.nc")


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

    

    


