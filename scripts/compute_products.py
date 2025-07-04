"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.
"""

# Ensure the top level directory has been added to PYTHONPATH
import argparse

#import functions 
import os
import yaml
from yaml.loader import SafeLoader


#import needed local functions
from osop.compute_products_func import calc_products



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
        help="variable to download, 2m_temperature, total_precipitation",
    )
    parser.add_argument("--downloaddir", required=True, help="location to get from")
    parser.add_argument("--productsdir", required=True, help="location to download to")
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
    Call the main funcitons to calculate anomalies and terciles
    """

    # get command line args
    args = parse_args()

    # unpack args and reformat if needed
    centre = args.centre
    downloaddir = args.downloaddir
    productsdir = args.productsdir
    month = int(args.month)
    leads = args.leads
    leadtime_month = [int(l) for l in args.leads.split(",")]
    leads_str = "".join([str(mon) for mon in leadtime_month])
    area = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    variable = args.variable

    # get remaning arguments from yml file
    ymllocation = os.path.join(downloaddir, "parseyml.yml")

    with open(ymllocation, "r") as stream:
        try:
            # Converts yaml document to python object
            services = yaml.load(stream, Loader=SafeLoader)
            # Converts contents to useable dictionary
            Services = services["Services"]
        except yaml.YAMLError as e:
            print(e)

    # add arguments to config
    config = dict(
        start_month=month,
        leads=leadtime_month,
        origin=centre,
        area_str=area_str,
        leads_str=leads_str,
        var=variable,
    )

    if args.years:
        config["hcstarty"] = args.years[0]
        config["hcendy"] = args.years[1]
    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016

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
