# Libraries for working with multi-dimensional arrays
import xarray as xr
import pandas as pd
import numpy as np
import os

# Forecast verification metrics with xarray
import xskillscore as xs

# Disable warnings for data download via API and matplotlib (do I need both???)
import warnings
warnings.filterwarnings('ignore')

from compute_products import get_thresh

import argparse
import sys

sys.path.append('/home/h01/edyer/IASAS/osop/')
from lib.osop.constants import SYSTEMS


def read_obs(obs_fname, config):
    era5_1deg = xr.open_dataset(obs_fname, engine='cfgrib')

    # Renaming to match hindcast names 
    era5_1deg = era5_1deg.rename({'latitude':'lat','longitude':'lon','time':'start_date'}).swap_dims({'start_date':'valid_time'})

    # Assign 'forecastMonth' coordinate values
    fcmonths = [mm+1 if mm>=0 else mm+13 for mm in [t.month - config['start_month'] for t in pd.to_datetime(era5_1deg.valid_time.values)] ]
    era5_1deg = era5_1deg.assign_coords(forecastMonth=('valid_time',fcmonths))
    # Drop obs values not needed (earlier than first start date) - this is useful to create well shaped 3-month aggregations from obs.
    era5_1deg = era5_1deg.where(era5_1deg.valid_time>=np.datetime64('{hcstarty}-{start_month:02d}-01'.format(**config)),drop=True)

    # CALCULATE 3-month AGGREGATIONS
    # NOTE rolling() assigns the label to the end of the N month period
    print('Calculate observation 3-monthly aggregations')
    # NOTE care should be taken with the data available in the "obs" xr.Dataset so the rolling mean (over valid_time) is meaningful
    era5_1deg_3m = era5_1deg.rolling(valid_time=3).mean()
    era5_1deg_3m = era5_1deg_3m.where(era5_1deg_3m.forecastMonth>=3)

    # As we don't need it anymore at this stage, we can safely remove 'forecastMonth'
    era5_1deg = era5_1deg.drop('forecastMonth')
    era5_1deg_3m = era5_1deg_3m.drop('forecastMonth')

    return era5_1deg, era5_1deg_3m


def scores_dtrmnstc(era5_1deg, era5_1deg_3m, hcst_bname, downloaddir):
    # Loop over aggregations
    for aggr in ['1m','3m']:

        if aggr=='1m':
            o = era5_1deg
        elif aggr=='3m':
            o = era5_1deg_3m
        else:
            raise BaseException(f'Unknown aggregation {aggr}')

        print(f'Computing deterministic scores for {aggr}-aggregation')

        # Read anomalies file
        h = xr.open_dataset(f'{downloaddir}/{hcst_bname}.{aggr}.anom.nc' )
        is_fullensemble = 'number' in h.dims

        l_corr=list()
        l_corr_pval=list()

        for this_fcmonth in h.forecastMonth.values:
            print(f'forecastMonth={this_fcmonth}' )
            thishcst = h.sel(forecastMonth=this_fcmonth).swap_dims({'start_date':'valid_time'})
            print(thishcst)
            thisobs = o.where(o.valid_time==thishcst.valid_time,drop=True)
            print(thisobs)
            thishcst_em = thishcst if not is_fullensemble else thishcst.mean('number')
            l_corr.append( xs.spearman_r(thishcst_em, thisobs, dim='valid_time') )
            l_corr_pval.append ( xs.spearman_r_p_value(thishcst_em, thisobs, dim='valid_time') )

        print(f'Concatenating (by fcmonth) correlation for {aggr}-aggregation')
        corr=xr.concat(l_corr,dim='forecastMonth')
        corr_pval=xr.concat(l_corr_pval,dim='forecastMonth')

        print(f'Saving to netCDF file correlation for {aggr}-aggregation')    
        corr.to_netcdf(f'{downloaddir}/scores/{hcst_bname}.{aggr}.corr.nc')
        corr_pval.to_netcdf(f'{downloaddir}/scores/{hcst_bname}.{aggr}.corr_pval.nc')


def scores_prblstc(era5_1deg, era5_1deg_3m, hcst_bname, config, downloaddir):
    quantiles = [1/3., 2/3.]
    numcategories = len(quantiles)+1
    # Loop over aggregations
    #for aggr in ['1m','3m']:
    aggr = config['aggr']
    if aggr=='1m':
        o = era5_1deg
    elif aggr=='3m':
        o = era5_1deg_3m
    else:
        raise BaseException(f'Unknown aggregation {aggr}')

    print(f'Computing probabilistic scores for {aggr}-aggregation')
    
    # READ hindcast probabilities file
    probs_hcst = xr.open_dataset(f'{downloaddir}/{hcst_bname}.{aggr}.tercile_probs.nc')
    #probs_hcst = probs_hcst[config['var']]

    ## met office forecast has unexpected valid times e.g. not just JJA
    #if 'ukmo' in hcst_bname:
    #    print('ukmo')
    #    #probs_hcst['start_date'] = probs_hcst.start_date.dt.strftime('%Y-%m')
    #    probs_hcst = probs_hcst.where(probs_hcst.start_date.dt.month==config['start_month'],drop=True)
    #print(probs_hcst.coords['valid_time'].values)

    l_roc=list()
    l_rps=list()
    l_rocss=list()
    l_bs=list()
    l_rel=list()
    o['valid_time'] = o.valid_time.dt.strftime('%Y-%m')
    print(probs_hcst.forecastMonth.values)

    for this_fcmonth in probs_hcst.forecastMonth.values:
        print(f'forecastMonth={this_fcmonth}')
        thishcst = probs_hcst.sel(forecastMonth=this_fcmonth).swap_dims({'start_date':'valid_time'})
        # CALCULATE probabilities from observations
        print('We need to calculate probabilities (tercile categories) from observations')
        l_probs_obs=list()
        thishcst['valid_time'] = thishcst.valid_time.dt.strftime('%Y-%m')
        thiso = o.where(o.valid_time==thishcst.valid_time,drop=True)

        for icat in range(numcategories):

            o_lo,o_hi = get_thresh(icat, quantiles, thiso, dims=['valid_time'])
            probo = 1. * np.logical_and(thiso>o_lo, thiso<=o_hi)
            if 'quantile' in probo:
                probo=probo.drop('quantile')
            l_probs_obs.append(probo.assign_coords({'category':icat}))

        thisobs = xr.concat(l_probs_obs, dim='category')

        print('Now we can calculate the probabilistic (tercile categories) scores')
        thisroc = xr.Dataset()
        thisrps = xr.Dataset()
        thisrocss = xr.Dataset()
        thisbs = xr.Dataset()
        thisrel = xr.Dataset()
        for var in thishcst.data_vars:
            if var == 'tprate':
                var_obs = 'tp'
            else:
                var_obs = var
            thisroc[var] = xs.roc(thisobs[var_obs],thishcst[var], dim='valid_time', bin_edges=np.linspace(0,1,101))

            thisrps[var] = xs.rps(thisobs[var_obs],thishcst[var], dim='valid_time', category_edges=None, input_distributions='p')

            thisrocss[var] = (thisroc[var] - 0.5) / (1. - 0.5)
            bscat = list()
            relcat = list()
            for cat in thisobs[var_obs].category:

                thisobscat = thisobs[var_obs].sel(category=cat)
                thishcstcat = thishcst[var].sel(category=cat)
                
                bscat.append(xs.brier_score(thisobscat, thishcstcat, dim='valid_time'))
                relcat.append(xs.reliability(thisobscat.astype(int), thishcstcat))

            thisbs[var] = xr.concat(bscat,dim='category')
            thisrel[var] = xr.concat(relcat,dim='category')

        l_roc.append(thisroc)
        l_rps.append(thisrps)
        l_rocss.append(thisrocss)
        l_bs.append(thisbs)
        l_rel.append(thisrel)

    print('concat roc')
    roc=xr.concat(l_roc,dim='forecastMonth')
    print('concat rps')
    rps=xr.concat(l_rps,dim='forecastMonth')
    print('concat rocss')
    rocss=xr.concat(l_rocss,dim='forecastMonth')
    print('concat bs')
    bs=xr.concat(l_bs,dim='forecastMonth')
    print('concat rel')
    rel=xr.concat(l_rel,dim='forecastMonth')

    print('writing to netcdf rps')
    rps.to_netcdf(f'{downloaddir}/scores/{hcst_bname}.{aggr}.rps.nc')
    print('writing to netcdf bs')
    bs.to_netcdf(f'{downloaddir}/scores/{hcst_bname}.{aggr}.bs.nc')
    print('writing to netcdf roc')
    roc.to_netcdf(f'{downloaddir}/scores/{hcst_bname}.{aggr}.roc.nc')
    print('writing to netcdf rocss')
    rocss.to_netcdf(f'{downloaddir}/scores/{hcst_bname}.{aggr}.rocss.nc')
    print('writing to netcdf rel')
    rel.to_netcdf(f'{downloaddir}/scores/{hcst_bname}.{aggr}.rel.nc')


def calc_scores(config, downloaddir):
    hcst_bname = '{origin}_{system}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}'.format(**config)

    ## obs info
    obs_fname = '{fpath}/era5_{hcstarty}-{hcendy}_monthly_{start_month}_{obs_str}_{area_str}.grib'.format(fpath=downloaddir,**config)

    ## read obs
    era5_1deg, era5_1deg_3m = read_obs(obs_fname, config)
    ## calc scores
    #scores_dtrmnstc(era5_1deg, era5_1deg_3m, hcst_bname, downloaddir)
    scores_prblstc(era5_1deg, era5_1deg_3m, hcst_bname, config, downloaddir)


def parse_args():
    """
    set up argparse to get command line arguments

    Returns:
        args: argparse args object
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--centre", required=True, help="centre to download")
    parser.add_argument("--month", required=True, help="start month for hindcasts")
    parser.add_argument("--aggregation", required=True, help="1m or 3m period for verification")
    parser.add_argument("--variable", required=True, help="variable to verify")
    parser.add_argument(
        "--leads", required=True, help="forecast range in months (comma separated)"
    )
    parser.add_argument(
        "--area",
        required=True,
        help="sub-area in degrees for retrieval (comma separated N,W,S,E)",
    )
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
    calculation of verification metrics
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
    obs_str = "".join([str(mon-1) for mon in leadtime_month])
    area = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    aggr = args.aggregation
    var = args.variable

    # add arguments to config
    config = dict(
        start_month = month,
        origin = centre,
        area_str = area_str,
        leads_str = leads_str,
        obs_str = obs_str,
        aggr = aggr,
        var = var
    )

    if args.years:
        config['hcstarty'] = args.years[0]
        config['hcendy'] = args.years[1]
    else:
        config['hcstarty'] = 1993
        config['hcendy'] = 2016

    # creat directory for verification scores
    scores_dir = os.path.join(downloaddir, 'scores')
    if not os.path.exists(scores_dir):
        os.makedirs(scores_dir)
            
    # hindcast info
    if centre == "eccc":
        # two models aka systems are live - call twice with each system number
        config['system'] = SYSTEMS["eccc_can"]
        calc_scores(config, downloaddir)
    
        ## repeat for second system
        config['system'] = SYSTEMS["eccc_gem5"]
        calc_scores(config, downloaddir)
    else:
        if centre not in SYSTEMS.keys():
            raise ValueError(f"Unknown system for C3S: {centre}")
        config['system'] = SYSTEMS[centre]
        calc_scores(config, downloaddir)