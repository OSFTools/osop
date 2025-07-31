"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.
"""

# Ensure the top level directory has been added to PYTHONPATH


#import functions 
import os
import yaml
from yaml.loader import SafeLoader


#import needed local functions
from osop.compute_products_func import calc_products


if __name__ == "__main__":
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the main funcitons to calculate anomalies and terciles
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
            productsdir = config_test["productsdir"]
            years = config_test["years"]
            centre = config_test["centre"]
            print('yml success')
        except yaml.YAMLError as e:
            print(e)
    
    downloaddir_test = os.path.expandvars(downloaddir_test)
    productsdir = os.path.expandvars(productsdir)
    os.makedirs(productsdir, exist_ok=True)


    # unpack args and reformat if needed
    downloaddir = downloaddir_test
    productsdir = productsdir
    month = int(Month_test)
    leads = leads_test
    leadtime_month = [int(l) for l in leads_test.split(",")]
    leads_str = "".join([str(mon) for mon in leadtime_month])
    area = [float(pt) for pt in area_test.split(",")]
    area_str = area_test.replace(",", ":")
    variable = varaible_test

   
    # add arguments to config
    config = dict(
        start_month=month,
        leads=leadtime_month,
        area_str=area_str,
        leads_str=leads_str,
        var=variable,
    )

    if years:
        config["hcstarty"] = years[0]
        config["hcendy"] = years[1]
    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016

    for centre in config_test["centre"]:
        config["origin"] = centre
        ## hindcast info
        if centre == "eccc":
            # two models aka systems are live - call twice with each system number
            config["system"] = Services["eccc_can"]
            calc_products(config, downloaddir, productsdir)

            ## repeat for second system
            config["system"] = Services["eccc_gem5"]
            calc_products(config, downloaddir, productsdir)
        else:
            if centre not in Services.keys():
                raise ValueError(f"Unknown system for C3S: {centre}")
            config["system"] = Services[centre]
            calc_products(config, downloaddir, productsdir)
