# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Script to compute scores for a given set of hindcasts and observations.

:note:
    This script currently can only be used to create scores using ERA5 as the comparison dataset
"""

# Ensure the top level directory has been added to PYTHONPATH
import argparse
from datetime import datetime
import logging
import os

import yaml
from yaml.loader import SafeLoader

# import local modules for function usage
from osop.compute_scores_func import calc_scores

logger = logging.getLogger(__name__)


def parse_args():
    """Set up argparse to get command line arguments.

    Returns
    -------
        args: argparse args object
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--centre", required=True, help="centre to download")
    parser.add_argument(
        "--obs_dataset",
        required=False,
        help="name of observation or reanalysis dataset. Optional. Defaults to ERA5",
    )
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
    parser.add_argument(
        "--downloaddir", required=True, help="location to get grib from"
    )
    parser.add_argument("--scoresdir", required=True, help="location to download to")
    parser.add_argument(
        "--productsdir", required=True, help="location to get products from"
    )
    parser.add_argument(
        "--years",
        required=False,
        help="Years to rerieve data for (comma separated). Optional. Default is hindcast period 1993-2016.",
    )
    parser.add_argument("--logdir", required=True, help="location to store logfiles")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the main function to do the actual
    calculation of verification metrics
    """

    # get command line args
    args = parse_args()

    # start logging - need to know logdir location before we can set it up
    logfile = os.path.join(
        args.logdir,
        f"scores_log_{args.variable}_{args.centre}_{args.month}_{datetime.today().strftime('%Y-%m-%d_%H:%M:%S')}.txt",
    )
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

    # unpack args and reformat if needed
    centre = args.centre
    downloaddir = args.downloaddir
    scoresdir = args.scoresdir
    productsdir = args.productsdir
    month = int(args.month)
    leads = args.leads
    leadtime_month = [int(l) for l in args.leads.split(",")]
    leads_str = "".join([str(mon) for mon in leadtime_month])
    obs_month = [int(l) - 1 for l in args.leads.split(",")]
    obs_str = "".join([str(mon) for mon in obs_month])
    area = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    hc_var = args.variable

    if hc_var == "2m_temperature":
        var = "t2m"
    elif hc_var == "total_precipitation":
        var = "tprate"
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
            logger.error(f"Error reading YAML file: {e}", stack_info=True)

    if args.years:
        config["hcstarty"] = args.years[0]
        config["hcendy"] = args.years[1]
    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016

    if args.obs_dataset:
        config["obs_name"] = args.obs
    else:
        config["obs_name"] = "era5"

    logger.debug(config)
    # hindcast info
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
