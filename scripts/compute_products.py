# Libraries for working with multi-dimensional arrays
import xarray as xr
import pandas as pd
import numpy as np
import os

# Date and calendar libraries
from dateutil.relativedelta import relativedelta

# Disable warnings for data download via API
import warnings
warnings.filterwarnings('ignore')

import argparse
import sys

sys.path.append('/home/h01/edyer/IASAS/osop/')
from lib.osop.constants import SYSTEMS


def calc_anoms(hcst_fname, hcst_bname, config, st_dim_name, DATADIR):
    print('Reading HCST data from file')
    hcst = xr.open_dataset(hcst_fname,engine='cfgrib', backend_kwargs=dict(time_dims=('forecastMonth', st_dim_name)))
    # force dask.array using chunks on leadtime, latitude and longitude coordinate
    hcst = hcst.chunk({'forecastMonth':1, 'latitude':'auto', 'longitude':'auto'})
    hcst = hcst.rename({'latitude':'lat','longitude':'lon', st_dim_name:'start_date'})

    print ('Re-arranging time metadata in xr.Dataset object')
    # Add start_month to the xr.Dataset
    start_month = pd.to_datetime(hcst.start_date.values[0]).month
    hcst = hcst.assign_coords({'start_month':start_month})
    # Add valid_time to the xr.Dataset
    vt = xr.DataArray(dims=('start_date','forecastMonth'), coords={'forecastMonth':hcst.forecastMonth,'start_date':hcst.start_date})
    vt.data = [[pd.to_datetime(std)+relativedelta(months=fcmonth-1) for fcmonth in vt.forecastMonth.values] for std in vt.start_date.values]
    hcst = hcst.assign_coords(valid_time=vt)

    # CALCULATE 3-month AGGREGATIONS
    # NOTE rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
    print('Computing 3-month aggregation')
    # rollng method defaults to look backwards
    hcst_3m = hcst.rolling(forecastMonth=3).mean()
    # Want only 3 month mean with complete 3 months
    hcst_3m = hcst_3m.where(hcst_3m.forecastMonth>=int(config['leads'][2]),drop=True)
    
    # CALCULATE ANOMALIES (and save to file)
    print('Computing anomalies 1m')
    hcmean = hcst.mean(['number','start_date'])
    anom = hcst - hcmean
    anom = anom.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))

    print('Computing anomalies 3m')
    hcmean_3m = hcst_3m.mean(['number','start_date'])
    anom_3m = hcst_3m - hcmean_3m
    anom_3m = anom_3m.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))

    print('Saving anomalies 1m/3m to netCDF files')
    anom.to_netcdf(f'{DATADIR}/{hcst_bname}.1m.anom.nc')
    anom_3m.to_netcdf(f'{DATADIR}/{hcst_bname}.3m.anom.nc')

    return hcst, hcst_3m


def get_thresh(icat,quantiles,xrds,dims=['number','start_date']):
    ''' Function to calculate the
    boundaries of forecast categories defined by quantiles'''

    if not all(elem in xrds.dims for elem in dims):           
        raise Exception('Some of the dimensions in {} is not present in the xr.Dataset {}'.format(dims,xrds)) 
    else:
        if icat == 0:
            xrds_lo = -np.inf
            xrds_hi = xrds.quantile(quantiles[icat],dim=dims)      
            
        elif icat == len(quantiles):
            xrds_lo = xrds.quantile(quantiles[icat-1],dim=dims)
            xrds_hi = np.inf
            
        else:
            xrds_lo = xrds.quantile(quantiles[icat-1],dim=dims)
            xrds_hi = xrds.quantile(quantiles[icat],dim=dims)
      
    return xrds_lo,xrds_hi


def prob_terc(hcst_bname, hcst, hcst_3m, DATADIR):
    '''CALCULATE PROBABILITIES for tercile categories
    by counting members within each category'''

    print('Computing probabilities (tercile categories)')
    quantiles = [1/3., 2/3.]
    numcategories = len(quantiles)+1

    for aggr,h in [("1m",hcst), ("3m",hcst_3m)]:
        if os.path.exist(f'{DATADIR}/{hcst_bname}.{aggr}.tercile_probs.nc'):
            print("(f'{DATADIR}/{hcst_bname}.{aggr}.tercile_probs.nc') exists")
        else:
            print(f'Computing tercile probabilities {aggr}')

            l_probs_hcst=list()
            for icat in range(numcategories):

                h_lo,h_hi = get_thresh(icat, quantiles, h)
                probh = np.logical_and(h>h_lo, h<=h_hi).sum('number')/float(h.dims['number'])

                # Instead of using the coordinate 'quantile' coming from the hindcast xr.Dataset
                # we will create a new coordinate called 'category'
                if 'quantile' in probh:
                    probh = probh.drop('quantile')
                l_probs_hcst.append(probh.assign_coords({'category':icat}))

            print(f'Concatenating {aggr} tercile probs categories')
            probs = xr.concat(l_probs_hcst,dim='category')                    
            print(f'Saving {aggr} tercile probs netCDF files')
            probs.to_netcdf(f'{DATADIR}/{hcst_bname}.{aggr}.tercile_probs.nc')


def calc_products(config, downloaddir):
        hcst_bname = '{origin}_{system}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{var}'.format(**config)
        hcst_fname = f'{downloaddir}/{hcst_bname}.grib'

        # For the re-shaping of time coordinates in xarray.Dataset we need to select the right one 
        #  -> burst mode ensembles (e.g. ECMWF SEAS5) use "time". This is the default option
        #  -> lagged start ensembles (e.g. MetOffice GloSea6) use "indexing_time" (see CDS documentation about nominal start date)
        st_dim_name = 'time' if not config.get('isLagged',False) else 'indexing_time'

        ## calc anoms
        hcst, hcst_3m = calc_anoms(hcst_fname, hcst_bname, config, st_dim_name, downloaddir)
        ## calc terc probs
        prob_terc(hcst_bname, hcst, hcst_3m, downloaddir)


def parse_args():
    """
    set up argparse to get command line arguments

    Returns:
        args: argparse args object
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--centre", required=True, help="centre to download")
    parser.add_argument("--month", required=True, help="start month for hindcasts")
    parser.add_argument(
        "--leads", required=True, help="forecast range in months (comma separated)"
    )
    parser.add_argument(
        "--area",
        required=True,
        help="sub-area in degrees for retrieval (comma separated N,W,S,E)",
    )
    parser.add_argument(
        "--variable",
        required=True,
        help="variable to download, 2m_temperature, total_precipitation")
    parser.add_argument("--downloaddir", required=True, help="location to download to")
    parser.add_argument(
        "--years",
        required=False,
        help="Years to rerieve data for (comma separated). Optional. Default is hindcast period 1993-2016.",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the main funciton to do the actual
    calculation of anomalies and terciles
    """

    # get command line args
    args = parse_args()

    # unpack args and reformat if needed
    centre = args.centre
    downloaddir = args.downloaddir
    month = int(args.month)
    leads = args.leads
    leadtime_month = [int(l) for l in args.leads.split(",")]
    leads_str = "".join([str(mon) for mon in leadtime_month])
    area = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    variable = args.variable

    # add arguments to config
    config = dict(
        start_month = month,
        leads = leadtime_month,
        origin = centre,
        area_str = area_str,
        leads_str = leads_str,
        var = variable
    )

    if args.years:
        config['hcstarty'] = args.years[0]
        config['hcendy'] = args.years[1]
    else:
        config['hcstarty'] = 1993
        config['hcendy'] = 2016

    # TODO replace with Nick's code for determining whether lagged from Grib
    if centre == 'ukmo':
        config['isLagged'] = True

    ## hindcast info
    if centre == "eccc":
        # two models aka systems are live - call twice with each system number
        config['system'] = SYSTEMS["eccc_can"]
        calc_products(config, downloaddir)
    
        ## repeat for second system
        config['system'] = SYSTEMS["eccc_gem5"]
        calc_products(config, downloaddir)
    else:
        if centre not in SYSTEMS.keys():
            raise ValueError(f"Unknown system for C3S: {centre}")
        config['system'] = SYSTEMS[centre]
        calc_products(config, downloaddir)