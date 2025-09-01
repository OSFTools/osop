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
from osop.ens_plotting import plot_forecasts, plot_mme


# Ensure the top level directory has been added to PYTHONPATH
import argparse




def parse_args():
    """
    set up argparse to get command line arguments

    Returns:
        args: argparse args object
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--location", required=True, help="location; type None for all")
    parser.add_argument("--centre", required=True, help="centre to download")
    parser.add_argument("--month", required=True, help="start month for hindcasts")
    parser.add_argument(
        "--variable", required=True, help="variable to verify. t2m or tprate"
    )
    parser.add_argument(
        "--leads", required=True, help="forecast range in months (comma separated)"
    )
    parser.add_argument(
        "--area",
        required=True,
        help="sub-area in degrees for retrieval (comma separated N,W,S,E)",
    )
    parser.add_argument("--downloaddir", required=True, help="location to get yml from")
    parser.add_argument("--productsfcdir", required=True, help="location to get from")
    parser.add_argument("--plotsdir", required=True, help="location to save too")
    parser.add_argument(
        "--yearsfc",
        required=False,
        help="Years to rerieve data for (comma separated). Optional. Default is hindcast period 1993-2016.",
    )
    parser.add_argument("--years", required=False, help="location to save too")

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
    location = args.location
    centre = args.centre
    downloaddir = args.downloaddir
    productsfcdir = args.productsfcdir
    plotsdir = args.plotsdir
    month = int(args.month)
    leads = args.leads
    print(leads)
    leadtime_month = [int(l) for l in args.leads.split(",")]
    print(leadtime_month)
    leads_str = "".join([str(mon) for mon in leadtime_month])
    print(leads_str)
    obs_month = [int(l) - 1 for l in args.leads.split(",")]
    obs_str = "".join([str(mon) for mon in obs_month])
    area = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    hc_var = args.variable
    i = list(map(int, leadtime_month))
    i = [x-2 for x in i]

    if hc_var == "2m_temperature":
        var = "2m_temperature"
    elif hc_var == "total_precipitation":
        var = "total_precipitation"
    else:
        raise ValueError(f"Unknown hindcast variable: {hc_var}")

    # add arguments to config
    config = dict(
        start_month=month,
        origin=centre,
        area_str=area_str,
        leads_str=leads_str,
        leads=leadtime_month,
        obs_str=obs_str,
        var=var,
        hc_var=hc_var,
        i = i,
        border = location
    )
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
    

    if args.yearsfc:
        years = [int(yr) for yr in args.yearsfc.split(",")]
        config["fcstarty"] = years[0]
        config["fcendy"] = years[0]
    else:
        config["fcstarty"] = 1993
        config["fcendy"] = 2016
    
    
    for x in i:
        config["i"]=x
        if centre == "eccc":
            # two models aka systems are live - call twice with each system number
            config["systemfc"] = Services["eccc_can"]
            plot_forecasts(productsfcdir,plotsdir,config)

            ## repeat for second system
            config["systemfc"] = Services["eccc_gem5"]
            plot_forecasts(productsfcdir,plotsdir,config)
        elif centre == "mme":
            config["systemfc"] = Services["mme"]
            plot_mme(Services,config,productsfcdir, plotsdir)
        else:
            if centre not in Services.keys():
               raise ValueError(f"Unknown system for C3S: {centre}")
            config["systemfc"] = Services[centre]
            plot_forecasts(productsfcdir,plotsdir,config)
