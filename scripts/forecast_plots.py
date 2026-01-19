"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.

This script currently can only be used to create scores using ERA5 as the comparison dataset

"""

from datetime import datetime
import os
import yaml
from yaml.loader import SafeLoader
import logging

# import local modules for function usage
from osop.ens_plotting import plot_forecasts


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

    parser.add_argument("--logdir", required=True, help="location to save logs")

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

    # unpack args and reformat if needed.

    logfile = os.path.join(args.logdir, 
        f"plots_log_{args.variable}_{args.centre}_{args.month}_{datetime.today().strftime('%Y-%m-%d_%H:%M:%S')}.txt")
    
    loglev = logging.INFO  # can be an argument later if needed
    logging.basicConfig(
        level=loglev,
        filename=logfile,
        encoding="utf-8",
        filemode="w",
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M",
    )

    location = args.location
    centre = args.centre
    downloaddir = args.downloaddir
    productsfcdir = args.productsfcdir
    plotsdir = args.plotsdir
    month = int(args.month)
    leads = args.leads
    logging.info(f'Plotting FC, Centre: {centre}, Month: {month},'\
                 f' Leads: {leads}, Location: {location}')

    leadtime_month = [int(l) for l in args.leads.split(",")]
    leads_str = "".join([str(mon) for mon in leadtime_month])
    obs_month = [int(l) - 1 for l in args.leads.split(",")]
    obs_str = "".join([str(mon) for mon in obs_month])
    area = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    hc_var = args.variable
    i = list(map(int, leadtime_month))
    i = [x - 2 for x in i]

    if hc_var == "2m_temperature":
        var = "2m_temperature"
    elif hc_var == "total_precipitation":
        var = "total_precipitation"
    else:
        logging.error(f"Unknown hindcast variable: {hc_var}")
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
        i=i,
        border=location,
    )
    # get remaining arguments from yml file
    ymllocation = os.path.join(downloaddir, "parseyml.yml")

    with open(ymllocation, "r") as stream:
        try:
            services_doc = yaml.load(stream, Loader=SafeLoader)
            ServicesRaw = services_doc["Services"]


            # Convert Services back the original dictionary (service -> value)
            # Remove Weights 
            Services = {
                svc: (val[0] if isinstance(val, (list, tuple)) else val)
                for svc, val in ServicesRaw.items()
            }
        except yaml.YAMLError as e:
            logging.error(f"Error reading YAML file: {e}", stack_info=True)

    if args.yearsfc:
        years = [int(yr) for yr in args.yearsfc.split(",")]
        config["fcstarty"] = years[0]
        config["fcendy"] = years[0]
    else:
        config["fcstarty"] = 1993
        config["fcendy"] = 2016

    logging.debug(config)
    for x in i:
        config["i"] = x
        if centre == "eccc":
            # two models aka systems are live - call twice with each system number
            config["systemfc"] = Services["eccc_can"]
            plot_forecasts(productsfcdir, plotsdir, config)

            ## repeat for second system
            config["systemfc"] = Services["eccc_gem5"]
            plot_forecasts(productsfcdir, plotsdir, config)
        else:
            if centre not in Services.keys():
                logging.error(f"Unknown system for C3S: {centre}")
                raise ValueError(f"Unknown system for C3S: {centre}")
            config["systemfc"] = Services[centre]
            plot_forecasts(productsfcdir, plotsdir, config)
    logging.info("Completed forecast plots successfully")