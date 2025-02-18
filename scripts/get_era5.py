# CDS API
import cdsapi

import os

# Disable warnings for data download via API
import warnings
warnings.filterwarnings('ignore')

import argparse


def get_obs(obs_fname, config):
    """
    calls cdsapi for the requested period and area
    retrieves monthly averaged reanalysis ERA5 data

    Args:
        obs_fname(str): fname to save to, and to check not already downloaded
        config dictionary containing necessary arguments for cdsapi

    """
    c = cdsapi.Client()
    if os.path.exists(obs_fname):
        print(f'File {obs_fname} already exists')
        return obs_fname

    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
       {
            'product_type': 'monthly_averaged_reanalysis',
            "variable": [config['var'],],
            # NOTE from observations we need to go one year beyond so we have available all the right valid dates
            # e.g. Nov.2016 start date forecast goes up to April 2017             
            'year': ['{}'.format(yy) for yy in range(config['hcstarty'],config['hcendy']+2)],
            'month': ['{:02d}'.format((config['start_month']+leadm)%12) if config['start_month']+leadm!=12 else '12' for leadm in config['leads_obs']],
            'time': '00:00',
            # We can ask CDS to interpolate ERA5 to the same grid used by C3S seasonal forecasts
            'grid': '1/1',
            'area': config['area'],
            'format': 'grib',
        },
        obs_fname)
    return obs_fname


def parse_args():
    """
    set up argparse to get command line arguments

    Returns:
        args: argparse args object
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--month", required=True, help="start month for observations")
    parser.add_argument(
        "--leads_obs", required=True, help="observation range in months (comma separated)"
    )
    parser.add_argument(
        "--area",
        required=True,
        help="sub-area in degrees for retrieval (comma separated N,W,S,E)",
    )
    parser.add_argument("--downloaddir", required=True, help="location to download to")
    parser.add_argument("--variable", required=True, help="variable to download")
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
    Call the main funciton to do download era5
    """

    # get command line args
    args = parse_args()

    # unpack args and reformat if needed
    downloaddir = args.downloaddir
    month = int(args.month)
    leads_obs = args.leads_obs
    leadtime_month = [int(l) for l in args.leads_obs.split(",")]
    # for filename to keep consistent with hindcast filenames
    leads_str = "".join([str(mon) for mon in leadtime_month])
    area_bounds = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    var = args.variable

    # add arguments to config
    config = dict(
        start_month = month,
        leads_obs = leadtime_month,
        area = area_bounds,
        area_str = area_str,
        leads_str = leads_str,
        var = var,
    )

    if args.years:
        config['hcstarty'] = args.years[0]
        config['hcendy'] = args.years[1]
    else:
        config['hcstarty'] = 1993
        config['hcendy'] = 2016

    ## obs info
    obs_fname = '{fpath}/era5_{var}_{hcstarty}-{hcendy}_monthly_{start_month}_{leads_str}_{area_str}.grib'.format(fpath=downloaddir,**config)
    print(obs_fname)
    get_obs(obs_fname, config)