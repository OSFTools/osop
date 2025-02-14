# Libraries for working with multi-dimensional arrays
import xarray as xr
import pandas as pd
import numpy as np
import os

# Libraries for plotting and geospatial data visualisation
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Disable warnings for data download via API and matplotlib (do I need both???)
import warnings
from cartopy import crs as ccrs
import calendar
warnings.filterwarnings('ignore')

import argparse
import sys

sys.path.append('/home/h01/edyer/IASAS/osop/')
from lib.osop.constants import SYSTEMS

CATNAMES=['lower tercile', 'middle tercile', 'upper tercile']

def prep_titles(config):

    # PREPARE STRINGS for TITLES
    tit_line1 = '{origin} {system}'.format(**config)
    tit_line2_base = f'Start month: {calendar.month_abbr[config["start_month"]].upper()}'

    if config['aggr']=='1m':
        validmonth = config['valid_month']#config['start_month'] + (config['valid_month']-1)
        validmonth = validmonth if validmonth<=12 else validmonth-12
        tit_line2 = tit_line2_base + f' - Valid month: {calendar.month_abbr[validmonth].upper()}'
    elif config['aggr'] == '3m':
        #validmonths = [vm if vm<=12 else vm-12 for vm in [config['start_month'] + (config['valid_month']-1) - shift for shift in range(3)]]
        validmonths = [vm if vm<=12 else vm-12 for vm in [config['valid_month'] + shift for shift in range(3)]]
        validmonths = [calendar.month_abbr[vm][0] for vm in validmonths] #reversed(validmonths)]
        tit_line2 = tit_line2_base + f' - Valid months: {"".join(validmonths)}'
    else:
        raise BaseException(f'Unexpected aggregation')
    return tit_line1, tit_line2


def plot_score(score_f, category, config, score, titles, datadir):

    # choose a projection and set up an axes based on that
    projection = ccrs.Mercator()
    ax = plt.axes(projection=projection)
    # Specify CRS
    crs = ccrs.PlateCarree()

    info = f'{config["var"]}_{config["origin"]}_{config["system"]}_fcst_start_month_{config["start_month"]}_leads_{config["leads_str"]}_{config["aggr"]}_cat_{category}_score_{score}'

    if score == 'rps':
        p = score_f[config['var']][0,:,:]
        lon = score_f[config['var']].lon
        lat = score_f[config['var']].lat
        cols = 'YlGn_r'
        ex_dir = 'max'
        under = 'purple'
        levels = np.linspace(0.,0.5,11)
        info = f'{config["var"]}_{config["origin"]}_{config["system"]}_fcst_start_month_{config["start_month"]}_leads_{config["leads_str"]}_{config["aggr"]}_score_{score}'
        plt.title(f'{score} \n' + titles[0] + f' {config["var"]}\n' + titles[1], loc='left')
    elif score == 'bs':
        p = score_f[config['var']].sel(category=category)[0,:,:]

        lon = score_f[config['var']].sel(category=category).lon
        lat = score_f[config['var']].sel(category=category).lat
        cols = 'YlGn_r'
        ex_dir = 'max'
        under = 'purple'
        plt.title(f'{score} \n' + titles[0] + f' {config["var"]}' + f' ({CATNAMES[category]})\n' + titles[1], loc='left')
        levels = np.linspace(0.,0.5,11)
    elif score in ['roc', 'rocss']:
        plt.title(f'{score} \n' + titles[0] + f' {config["var"]}' + f' ({CATNAMES[category]})\n' + titles[1], loc='left')
        p = score_f[config['var']].sel(category=category)[0,:,:]
        lon = score_f[config['var']].sel(category=category).lon
        lat = score_f[config['var']].sel(category=category).lat
        cols = 'YlGn'
        ex_dir = 'min'
        under = 'lightgray'
        levels = np.linspace(0.5,1.,6)

    cs = plt.contourf(lon, lat, p, transform=ccrs.PlateCarree(),levels=levels,cmap=cols, extend=ex_dir)
    cs.cmap.set_under(under)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2.)
   
    print(info)
    plt.colorbar()
    plt.savefig(os.path.join(datadir, 'scores', f'{info}.png'))
    plt.close()


def plot_rel(score_f, category, config, score, datadir):
    info = f'{config["var"]}_{config["origin"]}_{config["system"]}_fcst_start_month_{config["start_month"]}_{config["aggr"]}_cat_{category}_score_{score}'
    
    for var in score_f.data_vars:
        print(score_f)
        p = score_f[var][0,:,:]

        fig, ax1 = plt.subplots(figsize=(18,10))
        #ax2 = ax1.twinx()

        ax1.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Perfect Reliability')
        for icat in p.category.values:
            cat_lab = f' {CATNAMES[icat]}'
            avalues = p.sel(category=icat)
            ax1.plot(avalues.forecast_probability.values, avalues.values, label=cat_lab)
            #ax2.hist(avalues.values, bins=10, alpha=0.5, label=f'category {icat}')

        ax1.set_xlabel('Forecast Probability')
        ax1.set_ylabel('Values')
        #ax2.set_ylabel('Histogram')

        ax1.legend(loc='upper left')
        #ax2.legend(loc='upper right')

        plt.tight_layout()  
        plt.savefig(os.path.join(datadir, 'scores', f'{info}.png'))
        plt.close()


def corr_plots(datadir, hcst_bname, aggr, config, titles):
    # Read the data files
    corr = xr.open_dataset(f'{datadir}/scores/{hcst_bname}.{aggr}.corr.nc')
    corr_pval = xr.open_dataset(f'{datadir}/scores/{hcst_bname}.{aggr}.corr_pval.nc')
    # RE-ARRANGE the DATASETS longitude values for plotting purposes
    corr = corr.assign_coords(lon=(((corr.lon + 180) % 360) - 180)).sortby('lon')
    corr_pval = corr_pval.assign_coords(lon=(((corr_pval.lon + 180) % 360) - 180)).sortby('lon')


    fig = plt.figure(figsize=(18,10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2.)
    corrvalues = corr[config['var']][0,:,:].values 
    corrpvalvalues = corr_pval[config['var']][0,:,:].values

    if corrvalues.T.shape == (corr[config['var']].lat.size, corr[config['var']].lon.size):
        print('Data values matrices need to be transposed')
        corrvalues = corrvalues.T
        corrpvalvalues = corrpvalvalues.T
    elif corrvalues.shape == (corr[config['var']].lat.size, corr[config['var']].lon.size):
        pass                           
        # print('Data values matrices shapes are ok')
    else:
        raise BaseException(f'Unexpected data value matrix shape: {corrvalues.shape}' )
    plt.contourf(corr[config['var']].lon,corr[config['var']].lat,corrvalues,levels=np.linspace(-1.,1.,11),cmap='RdYlBu_r')
    cb = plt.colorbar(shrink=0.5)
    cb.ax.set_ylabel('Correlation',fontsize=12)
    origylim = ax.get_ylim()
    plt.contourf(corr_pval[config['var']].lon,corr_pval[config['var']].lat,corrpvalvalues,levels=[0.05,np.inf],hatches=['...',None],colors='none')
    # We need to ensure after running plt.contourf() the ylim hasn't changed
    if ax.get_ylim()!=origylim:
        ax.set_ylim(origylim)

    plt.title(f'{titles[0]} (stippling where significance below 95%)\n {titles[1]}',loc='left')
    plt.tight_layout()
    figname = f'{datadir}/scores/{hcst_bname}.{aggr}.corr.png'
    plt.savefig(figname)

def generate_plots(config, titles, downloaddir):
    ## read in the data
    
    score_fname = '{origin}_{system}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{fname_var}.{aggr}.{score}.nc'.format(**config)
    score_data = xr.open_dataset(os.path.join(downloaddir, 'scores', score_fname))
    if config['score'] == 'corr':
        score_fname = '{origin}_{system}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{fname_var}'.format(**config)
        corr_plots(downloaddir, score_fname, config['aggr'], config, titles)

    elif config['score'] == 'rel':
        plot_rel(score_data, None, config, config['score'], downloaddir)
    
    elif config['score'] == 'rps':
        plot_score(score_data, None, config, config['score'], titles, downloaddir)
    
    elif config['score'] == 'bs':
        for cat in [0,1,2]:
            plot_score(score_data, cat, config, config['score'], titles, downloaddir)

    else:
        for cat in [0,1,2]:
            config['category'] = cat
            plot_score(score_data, cat, config, config['score'], titles, downloaddir)

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
    parser.add_argument(
        "--variable",
        required=True,
        help="variable to verify. 2m_temperature or total_precipitation")
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


if __name__ == '__main__':
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the plotting functions to generate verification plots
    """
    scores = ['corr', 'roc', 'rocss', 'rps', 'rel', 'bs']

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
    fname_var = args.variable

    if fname_var == "2m_temperature":
        var = "t2m"
    elif fname_var == "total_precipitation":
        var = "tprate"
    else:
        raise ValueError(f"Unknown hindcast variable: {fname_var}")

    valid_month = month + (leadtime_month[0]-1)

    # add arguments to config
    config = dict(
        start_month = month,
        valid_month = valid_month,
        origin = centre,
        area_str = area_str,
        leads_str = leads_str,
        obs_str = obs_str,
        aggr = aggr,
        fname_var = fname_var,
        var = var
    )

    if args.years:
        config['hcstarty'] = args.years[0]
        config['hcendy'] = args.years[1]
    else:
        config['hcstarty'] = 1993
        config['hcendy'] = 2016

 
    for score in scores:
        config['score'] = score

        # run for appropriate system
        if centre == "eccc":
            # two models aka systems are live - call twice with each system number
            config['system'] = SYSTEMS["eccc_can"]
            ## set titles
            titles = prep_titles(config)
            generate_plots(config, titles, downloaddir)
        
            ## repeat for second system
            config['system'] = SYSTEMS["eccc_gem5"]
            ## set titles
            titles = prep_titles(config)
            generate_plots(config, titles, downloaddir)
        else:
            if centre not in SYSTEMS.keys():
                raise ValueError(f"Unknown system for C3S: {centre}")
            config['system'] = SYSTEMS[centre]
            ## set titles
            titles = prep_titles(config)
            generate_plots(config, titles, downloaddir)