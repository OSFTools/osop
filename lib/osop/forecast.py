import xarray as xr
import xarray as xr
import pandas as pd
import eccodes
import matplotlib.pyplot as plt

## Import The files for practice: the data is using 
# Month: 5 - May
# Leads: 2,3,4  - June, July, August
# Variable: t2m 
# Services: 9

print("START")
#path to the .nc file
MF_HC_TT = '/data/scratch/eleanor.dean/seafoam/data/master/hindcast/products/meteo_france_9_1993-2016_monthly_mean_5_234_45:-30:-2.5:60_total_precipitation.3m.tercile_thresholds.nc'
MF_HC_ = '/data/scratch/eleanor.dean/seafoam/data/master/hindcast/downloads/meteo_france_9_1993-2016_monthly_mean_5_234_45:-30:-2.5:60_total_precipitation.grib'
MF_FC_ = '/data/scratch/eleanor.dean/seafoam/data/master/forecast/downloads/meteo_france_9_2025-2025_monthly_mean_5_234_45:-30:-2.5:60_total_precipitation.grib'
# Open the file as xr dataset
MF_HC_TT = xr.open_dataset(MF_HC_TT)
MF_HC = xr.open_dataset(MF_HC_)
MF_FC = xr.open_dataset(MF_FC_)
#Check the dataset info is correct
print("this is hctt:",MF_HC_TT)
# print("this is hc:",MF_HC)
# print("this is fc:",MF_FC)

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

## Check all the data before working on it
print("this is FC", FC)
print(FC.forecastMonth.values)

"""
Okay so the my plan is to: select a month - June, then Select each catagory of the TC array for this month. 
Then for each EM create a boolean mask for which catagory it falls into. Concat that back into one array. Sum along the Number for True. 
Turn that value into a percentage. Plot a small map for check -- relevant to 1m code
"""
###########################################################So this is for 3m code #######################################
#print("this is FC:",FC.tprate.values)
fcst_3m = FC.rolling(forecastMonth=3).mean()
#print("this is fcst_3m", fcst_3m)
# # Want only 3 month mean with complete 3 months
fcst_3m = fcst_3m.isel(forecastMonth=(fcst_3m['forecastMonth'] == 4))

# print("this is fcst_3m with drop", fcst_3m)
# print("this is FC", FC)
# print("this is hctt,", MF_HC_TT)
data_var = list(FC.data_vars)[0]
tprate = FC[data_var]
#FC_SUM = FC.mean(dim='forecastMonth')
#FC_SUM = FC.rolling(forecastMonth=3).mean()
#FC_SUM = FC.sum(dim='forecastMonth' )
#FC_SUM = FC_SUM.drop(labels='surface')
#FC_SUM = FC.sum(dim='forecastMonth' )
# FC_SUMM = list_test.append(FC_SUM)
print("this is test")
#print(FC_SUM)
#Finsihed at the case where this seems to work for plan but needs to return to an x array that is getting added with each section



######### This is where the edit will be for total precip
fc=fcst_3m
hctt = MF_HC_TT

#hcmean_3m = hcst_3m.mean(["number", "start_date"])
#fc_ens_mean_3m = fcst_3m.mean(["number"])
##anom_3m = hcst_3m - hcmean_3m
##anom_3m = anom_3m.assign_attrs(reference_period="{hcstarty}-{hcendy}".format(**config))

data_var = list(fc.data_vars)[0]
t2mfc = fc[data_var]
print("this is t2mfc", t2mfc.values)
data_var = list(hctt.data_vars)[0]
print("this is size for one ensemble member:", t2mfc.isel(number = 1).values.size)
print("this is hctt.t2m size:",hctt.tprate.values)

### Select catagory
cat_0_thresh = hctt[data_var].sel(category=0)
print("this is the size for one catagory:", cat_0_thresh.values.size)
cat_1_thresh = hctt[data_var].sel(category=1)

#file = '/data/scratch/eleanor.dean/seafoam/data/master/forecast/products/meteo_france_9_2025-2025_monthly_mean_5_234_45:-30:-2.5:60_2m_temperature.imonth_0.forecast_percentages.nc'

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
print(true_counts.values)
true_counts1 = cat_1_mask_array.sum(dim='number')
true_counts3 = cat_3_mask_array.sum(dim='number')
print("this is true_counts;", true_counts, true_counts.values)

### Calculate that as a percentage 
percentage_true = (true_counts / cat_0_mask_array.sizes['number']) * 100
percentage_true1 = (true_counts1 / cat_1_mask_array.sizes['number']) * 100
percentage_true3 = (true_counts3 / cat_3_mask_array.sizes['number']) * 100

 ### Sanity check
print(percentage_true.values)

#compare the dataset from compare terciles 
file = '/data/scratch/eleanor.dean/seafoam/data/master/forecast/products/meteo_france_9_2025-2025_monthly_mean_5_234_45:-30:-2.5:60_total_precipitation.3m.forecast_percentages.nc'
compare = xr.open_dataset(file)
print(compare)
cat_0 = compare['lower']
print("this is the one that gets plotted:", percentage_true.values)
print("this is the import:", cat_0.values)

# ### Plot for second Sanity check

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import colors as c

data_to_plot = percentage_true1.values[0]
print("this is data to plot", data_to_plot)

lats= percentage_true1['lat'].values
lons = percentage_true1['lon'].values


fig = plt.figure(figsize=(20,15))
ax = plt.axes(projection=ccrs.PlateCarree())

ax.coastlines(resolution='50m')
ax.add_feature(cfeature.BORDERS, linestyle=':')

 #gl = ax.gridlines(drawlables=True, dms=True, x_inline=False, y_inline=False)
 #gl.top_labels = False
 #gl.right_labels = False
cMap = c.ListedColormap(['b','c','w','yellow','orange','darkorange', 'r'])

bounds = [0, 10, 20, 40, 50, 60, 70, 100]
norm = c.BoundaryNorm(bounds, cMap.N)

c = ax.pcolormesh(lons,lats, data_to_plot,cmap = cMap ,norm=norm, shading='auto')
cb = plt.colorbar(c, ax=ax, orientation='vertical', shrink=0.5, pad=0.05)

plt.show()



###########################################################this is where 3m code ends#######################################















###########################################################this is where 1m code starts#######################################
##Note if using this - make sure to change the import at the top to 1m 
# atm this is plotting upper tercile. '# check what month before use


# ###select June for both datasets#

# fc = FC.isel(forecastMonth=2)
# hctt = MF_HC_TT.isel(forecastMonth=2)
# #meteo_france_9_2025-2025_monthly_mean_5_234_45:-30:-2.5:60_2m_temperature.forecast_percentages.nc
# # Check the size for t2m for one ensemble meber agains catagory 
# data_var = list(fc.data_vars)[0]
# t2mfc = fc[data_var]

# data_var = list(hctt.data_vars)[0]
# print("this is size for one ensemble member:", t2mfc.isel(number = 1).values.size)
# print("this is hctt.t2m size:",hctt.t2m.size)

# ### Select catagory
# cat_0_thresh = hctt[data_var].sel(category=0)
# print("this is the size for one catagory:", cat_0_thresh.values.size)
# cat_1_thresh = hctt[data_var].sel(category=1)

# file = '/data/scratch/eleanor.dean/seafoam/data/master/forecast/products/meteo_france_9_2025-2025_monthly_mean_5_234_45:-30:-2.5:60_2m_temperature.imonth_0.forecast_percentages.nc'

# ### Create empty arrays for the mask to be added to
# cat_0_masks = []
# cat_1_masks = []
# cat_3_masks = []

# ### Go through each ensemble member (or number) and compare to catagory
# for i in range(t2mfc.sizes['number']):
#     t2m_slice = t2mfc.isel(number=i)
#     mask = t2m_slice < cat_0_thresh
#     cat_0_masks.append(mask)
#     mask1 = t2m_slice > cat_1_thresh
#     cat_1_masks.append(mask1)
#     mask3 = (t2m_slice > cat_0_thresh) & (t2m_slice < cat_1_thresh)
#     cat_3_masks.append(mask3)

# ### Concat each mask into a single DataArray via the ensemble members
# cat_0_mask_array = xr.concat(cat_0_masks, dim='number')
# cat_1_mask_array = xr.concat(cat_1_masks, dim='number')
# cat_3_mask_array = xr.concat(cat_3_masks, dim='number')
# # Check that thats worked appropriately - should be a datta array with coords
# # lon, lat, number and then be trues and falses 
# print("this is cat_0_mask_array:", cat_0_mask_array.values)

# ### Add up the trues along number for each catagory
# true_counts = cat_0_mask_array.sum(dim='number')
# true_counts1 = cat_1_mask_array.sum(dim='number')
# true_counts3 = cat_3_mask_array.sum(dim='number')
# print("this is true_counts;", true_counts, true_counts.values)

# ### Calculate that as a percentage 
# percentage_true = (true_counts / cat_0_mask_array.sizes['number']) * 100
# percentage_true1 = (true_counts1 / cat_1_mask_array.sizes['number']) * 100
# percentage_true3 = (true_counts3 / cat_3_mask_array.sizes['number']) * 100

#  ### Sanity check
# print(percentage_true.values)

# #compare the dataset from compare terciles 
# # file = '/data/scratch/eleanor.dean/seafoam/data/master/forecast/products/meteo_france_9_2025-2025_monthly_mean_5_234_45:-30:-2.5:60_2m_temperature.imonth_0.forecast_percentages.nc'
# # compare = xr.open_dataset(file)
# # print(compare)
# # cat_0 = compare['lower']
# # print("this is the one that gets plotted:", percentage_true.values)
# # print("this is the import:", cat_0.values)

# # ### Plot for second Sanity check

# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from matplotlib import colors as c

# data_to_plot = percentage_true1.values

# lats= percentage_true1['lat'].values
# lons = percentage_true1['lon'].values


# fig = plt.figure(figsize=(20,15))
# ax = plt.axes(projection=ccrs.PlateCarree())

# ax.coastlines(resolution='50m')
# ax.add_feature(cfeature.BORDERS, linestyle=':')

#  #gl = ax.gridlines(drawlables=True, dms=True, x_inline=False, y_inline=False)
#  #gl.top_labels = False
#  #gl.right_labels = False
# cMap = c.ListedColormap(['b','g','w','m','y','r'])

# c = ax.pcolormesh(lons,lats, data_to_plot,cmap = 'coolwarm' , shading='auto')
# cb = plt.colorbar(c, ax=ax, orientation='vertical', shrink=0.5, pad=0.05)

# plt.show()

###########################################################this is where 1m code ends#######################################

# ### Old code Snippets
# #cat_0_mask = t2mfc < cat_0_thresh
# #cat_1_mask = t2mfc > cat_1_thresh
# #cat_3_mask = (t2mfc > cat_0_thresh) & (t2mfc < cat_1_thresh)

# #cat_0_count = cat_0_mask.sum(dim='number')
# #cat_1_count = cat_1_mask.sum(dim='number')

# #total_members = t2mfc.sizes["number"]

# #cat_0_percent = (cat_0_count / total_members) * 100
# #cat_1_percent = (cat_1_count / total_members) * 100

# #print(cat_0_percent.values.size)
# #print(cat_1_percent.values.size)

##old code snippets for comparing the ones from this python file and the ones with the shell. 

# # results = "/data/scratch/eleanor.dean/seafoam/data/master/forecast/products/meteo_france_9_2025-2025_monthly_mean_5_234_45:-30:-2.5:60_2m_temperature.imonth_0.forecast_percentages.nc"
# # results = xr.open_dataset(results)


# # # Stack the three layers into a new dimension C
# # temperature_data = np.stack([
# #     results['lower'].values,
# #     results['middle'].values,
# #     results['higher'].values
# # ], axis=0)

# # # Define new coordinates
# # C = [1, 2, 3]  # 1: lower, 2: middle, 3: higher
# # Y = results['lat'].values
# # X = results['lon'].values

# # # Create the new dataset
# # new_dataset = xr.Dataset(
# #     {
# #         "temperature": (("C", "Y", "X"), temperature_data)
# #     },
# #     coords={
# #         "C": C,
# #         "Y": Y,
# #         "X": X
# #     }
# # )


# # #print(new_dataset)

# # #print(results)
# # #print(results['lower'].values)

# # from osop.ens_plotting import plot_tercile_fc
# # title = "meteofrance plot test"

# # output = plot_tercile_fc(new_dataset, title, l_borders=True, var="temperature", mask=None)

 
# # output.savefig('/data/scratch/eleanor.dean/seafoam/data/master/forecast/plots/my_plot.png')

