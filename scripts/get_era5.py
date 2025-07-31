"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.
"""

# CDS API
import cdsapi

import os

# Disable warnings for data download via API
import warnings
warnings.filterwarnings('ignore')


import yaml
from yaml.loader import SafeLoader


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


if __name__ == "__main__":
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the main funciton to do download era5
    """
    # get remaning arguments from yml file
    ymllocation = os.path.join("variables.yml")

    with open(ymllocation, "r") as stream:
        try:
            print('yml found')
            # Converts yaml document to python object
            config_test = yaml.load(stream, Loader=SafeLoader)
            # Converts contents to useable dictionary
            Services = config_test["Services"]
            Month_test = config_test["month"]
            leads_test = config_test["leads"]
            area_test = config_test["area"]
            location_test = config_test["location"]
            varaible_test = config_test["variable"]
            downloaddir_test = config_test["downloaddir"]
            years = config_test["years"]
            print('yml success')
        except yaml.YAMLError as e:
            print(e)

    #Make directory if it dosnt exist. 

    downloaddir_test = os.path.expandvars(downloaddir_test)
    os.makedirs(downloaddir_test, exist_ok=True)

    # unpack args and reformat if needed
    downloaddir = downloaddir_test
    month = Month_test
    leadtime_month = [int(l)-1 for l in leads_test.split(",")]
    # for filename to keep consistent with hindcast filenames
    leads_str = "".join([str(mon) for mon in leadtime_month])
    area_bounds = [float(pt) for pt in area_test.split(",")]
    area_str = area_test.replace(",", ":")
    var = varaible_test

    # add arguments to config
    config = dict(
        start_month = month,
        leads_obs = leadtime_month,
        area = area_bounds,
        area_str = area_str,
        leads_str = leads_str,
        var = var,
    )

    if years:
        config['hcstarty'] = int(years[0])
        config['hcendy'] = int(years[1])
    else:
        config['hcstarty'] = 1993
        config['hcendy'] = 2016

    ## obs info
    obs_fname = '{fpath}/era5_{var}_{hcstarty}-{hcendy}_monthly_{start_month}_{leads_str}_{area_str}.grib'.format(fpath=downloaddir,**config)
    print(obs_fname)
    get_obs(obs_fname, config)