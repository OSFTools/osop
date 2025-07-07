import xarray as xr
import xarray as xr
import pandas as pd
import numpy as np
import eccodes
import matplotlib.pyplot as plt




### Import The files for practice: the data is using 
# Month: 5 - May
# Leads: 2,3,4  - June, July, August
# Variable: t2m 
# Services: 9

#path to the .nc file
MF_HC_TT = '/data/scratch/eleanor.dean/seafoam/data/master/hindcast/products/meteo_france_9_1993-2016_monthly_mean_5_234_45:-30:-2.5:60_2m_temperature.1m.tercile_thresholds.nc'
MF_HC_ = '/data/scratch/eleanor.dean/seafoam/data/master/hindcast/downloads/meteo_france_9_1993-2016_monthly_mean_5_234_45:-30:-2.5:60_2m_temperature.grib'
MF_FC_ = '/data/scratch/eleanor.dean/seafoam/data/master/forecast/downloads/meteo_france_9_2025-2025_monthly_mean_5_234_45:-30:-2.5:60_2m_temperature.grib'
# Open the file as xr dataset
MF_HC_TT = xr.open_dataset(MF_HC_TT)
MF_HC = xr.open_dataset(MF_HC_)
MF_FC = xr.open_dataset(MF_FC_)
#Check the dataset info is correct
print("this is hctt:",MF_HC_TT.t2m.size)
print("this is hc:",MF_HC)
print("this is fc:",MF_FC)


### STOLEN FUNCTIONS from compute_products -- reformates the grib FC to be the same as the TC

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

st_dim_name = get_tindex(MF_FC_)

def calc_anoms(hcst_fname, st_dim_name):
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

    print("Reading HCST data from file")
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

    print("Re-arranging time metadata in xr.Dataset object")
    # Add start_month to the xr.Dataset
    start_month = pd.to_datetime(hcst.start_date.values).month
    hcst = hcst.assign_coords({"start_month": start_month})
    # Add valid_time to the xr.Dataset
    return(hcst)

FC = calc_anoms(MF_FC_, st_dim_name)


### Check all the data before working on it
print("this is FC", FC)
print(FC.forecastMonth.values)


"""
Okay so the my plan is to: select a month - June, then Select each catagory of the TC array for this month. T
Then for each EM create a boolean mask for which catagory it falls into. Concat that back into one array. Sum along the Number for True. 
Turn that value into a percentage. Plot a small map for check
"""

###select June for both datasets#

fc = FC.isel(forecastMonth=0)
hctt = MF_HC_TT.isel(forecastMonth=0)

# Check the size for t2m for one ensemble meber agains catagory 
t2mfc = fc['t2m']
print("this is size for one ensemble member:", t2mfc.isel(number = 1).values.size)
print("this is hctt.t2m size:",hctt.t2m.size)


### Select catagory
cat_0_thresh = hctt['t2m'].sel(category=0)
print("this is the size for one catagory:", cat_0_thresh.values.size)
cat_1_thresh = hctt['t2m'].sel(category=1)


### Create empty arrays for the mask to be added to
cat_0_masks = []
cat_1_masks = []
cat_3_masks = []

### Go through each ensemble member (or number) and compare to catagory
for i in range(t2mfc.sizes['number']):
    t2m_slice = t2mfc.isel(number=i)
    mask = t2m_slice < cat_0_thresh
    cat_0_masks.append(mask)
    mask1 = t2m_slice > cat_1_thresh
    cat_1_masks.append(mask1)
    mask3 = (t2m_slice > cat_0_thresh) & (t2m_slice < cat_1_thresh)
    cat_3_masks.append(mask3)

### Concat each mask into a single DataArray via the ensemble members
cat_0_mask_array = xr.concat(cat_0_masks, dim='number')
cat_1_mask_array = xr.concat(cat_1_masks, dim='number')
cat_3_mask_array = xr.concat(cat_3_masks, dim='number')
# Check that thats worked appropriately - should be a datta array with coords
# lon, lat, number and then be trues and falses 
print("this is cat_0_mask_array:", cat_0_mask_array.values)

### Add up the trues along number for each catagory
true_counts = cat_0_mask_array.sum(dim='number')
true_counts1 = cat_1_mask_array.sum(dim='number')
true_counts3 = cat_3_mask_array.sum(dim='number')
print("this is true_counts;", true_counts, true_counts.values)

### Calculate that as a percentage 
percentage_true = (true_counts / cat_0_mask_array.sizes['number']) * 100
percentage_true1 = (true_counts1 / cat_1_mask_array.sizes['number']) * 100
percentage_true3 = (true_counts3 / cat_3_mask_array.sizes['number']) * 100


### Sanity check
print(percentage_true.values)
print(percentage_true1.values)


### Plot for second Sanity check

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import colors as c

data_to_plot = percentage_true.values

lats= percentage_true3['lat'].values
lons = percentage_true3['lon'].values



fig = plt.figure(figsize=(20,15))
ax = plt.axes(projection=ccrs.PlateCarree())

ax.coastlines(resolution='50m')
ax.add_feature(cfeature.BORDERS, linestyle=':')

#gl = ax.gridlines(drawlables=True, dms=True, x_inline=False, y_inline=False)
#gl.top_labels = False
#gl.right_labels = False
cMap = c.ListedColormap(['b','g','w','m','y','r'])



c = ax.pcolormesh(lons,lats, data_to_plot,cmap = 'coolwarm' , shading='auto')
cb = plt.colorbar(c, ax=ax, orientation='vertical', shrink=0.5, pad=0.05)


plt.show()


### Old code Snippets
#cat_0_mask = t2mfc < cat_0_thresh
#cat_1_mask = t2mfc > cat_1_thresh
#cat_3_mask = (t2mfc > cat_0_thresh) & (t2mfc < cat_1_thresh)

#cat_0_count = cat_0_mask.sum(dim='number')
#cat_1_count = cat_1_mask.sum(dim='number')

#total_members = t2mfc.sizes["number"]

#cat_0_percent = (cat_0_count / total_members) * 100
#cat_1_percent = (cat_1_count / total_members) * 100

#print(cat_0_percent.values.size)
#print(cat_1_percent.values.size)