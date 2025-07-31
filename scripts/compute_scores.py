"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.

This script currently can only be used to create scores using ERA5 as the comparison dataset

"""

import os
import yaml
from yaml.loader import SafeLoader


#import local modules for function usage 
from osop.compute_scores_func import calc_scores


# Ensure the top level directory has been added to PYTHONPATH


if __name__ == "__main__":

    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the main funciton to do the actual
    calculation of verification metrics
    """


    # get command line args
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
            varaible_test = config_test["variable"]
            downloaddir_test = config_test["downloaddir"]
            scoresdir = config_test["scoresdir"]
            productsdir = config_test["productsdir"]
            years = config_test["years"]
            centre = config_test["centre"]
            obs = config_test["obs"]
            print('yml success')
        except yaml.YAMLError as e:
            print(e)

    downloaddir_test = os.path.expandvars(downloaddir_test)
    scoresdir= os.path.expandvars(scoresdir)
    productsdir = os.path.expandvars(productsdir)
    os.makedirs(scoresdir, exist_ok=True)

    # unpack args and reformat if needed
    centre = centre
    downloaddir = downloaddir_test
    scoresdir = scoresdir
    productsdir = productsdir
    month = int(Month_test)
    leads = leads_test
    leadtime_month = [int(l) for l in leads_test.split(",")]
    leads_str = "".join([str(mon) for mon in leadtime_month])
    obs_month = [int(l) - 1 for l in leads_test.split(",")]
    obs_str = "".join([str(mon) for mon in obs_month])
    area = [float(pt) for pt in area_test.split(",")]
    area_str = area_test.replace(",", ":")
    hc_var = varaible_test

    if hc_var == "2m_temperature":
        var = "t2m"
    elif hc_var == "total_precipitation":
        var = "tprate"
    else:
        raise ValueError(f"Unknown hindcast variable: {hc_var}")

    # add arguments to config
    config = dict(
        start_month=month,
        area_str=area_str,
        leads_str=leads_str,
        leads=leadtime_month,
        obs_str=obs_str,
        var=var,
        hc_var=hc_var,
    )
    
    if years:
        config["hcstarty"] = years[0]
        config["hcendy"] = years[1]
    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016
    
    if obs:
        config['obs_name'] = obs
    else:
        config['obs_name'] = 'era5'


    # hindcast info
    for centre in config_test["centre"]:
        config["origin"] = centre
        if centre == "eccc":
            # two models aka systems are live - call twice with each system number
            config["system"] = Services["eccc_can"]
            calc_scores(config, downloaddir, scoresdir, productsdir)

            ## repeat for second system
            config["system"] = Services["eccc_gem5"]
            calc_scores(config, downloaddir, scoresdir, productsdir)
        else:
            if centre not in Services.keys():
                raise ValueError(f"Unknown system for C3S: {centre}")
            config["system"] = Services[centre]
            calc_scores(config, downloaddir, scoresdir, productsdir)
