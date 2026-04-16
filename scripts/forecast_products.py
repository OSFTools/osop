# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Script to calculate forecast products for a given centre, variable, month and lead time.

Notes
-----
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
from osop.compare_terciles import compute_forecast, mme_products

logger = logging.getLogger(__name__)


def parse_args():
    """Set up argparse to get command line arguments.

    Returns
    -------
        args: argparse args object
    """
    parser = argparse.ArgumentParser()
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
    parser.add_argument(
        "--downloaddir", required=True, help="location to get grib from"
    )
    parser.add_argument(
        "--downloadhcdir", required=True, help="location to get yml hc services from"
    )
    parser.add_argument(
        "--productshcdir", required=True, help="location to get products from"
    )
    parser.add_argument("--productsfcdir", required=True, help="location to save too")
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
    Call the main function to do the actual
    calculation of verification metrics
    """

    # get command line args
    args = parse_args()

    logfile = os.path.join(
        args.logdir,
        f"products_log_{args.variable}_{args.centre}_{args.month}_{datetime.today().strftime('%Y-%m-%d_%H:%M:%S')}.txt",
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
    downloadhcdir = args.downloadhcdir
    productshcdir = args.productshcdir
    productsfcdir = args.productsfcdir
    month = int(args.month)
    leads = args.leads
    logger.info(f"Doing FC products, Centre: {centre}, Month: {month}, Leads: {leads}")

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
        logger.error(f"Unknown hindcast variable: {hc_var}")
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

    def normalize_services(services_raw):
        """Services_raw(dict): config file from the passed through yml.

        Normalize a mapping of service -> value into service -> (id, weight)
        Accepts:
        - scalar (e.g., 51)                 -> (51, 1)
        - [id] or (id,)                     -> (id, 1)
        - [id, weight] or (id, weight, ...) -> (id, weight) (first two only)
        - []                                -> (None, 1) with a warning
        """
        if not isinstance(services_raw, dict):
            raise TypeError("'Services' must be a mapping of service -> [id, weight]")

        normalized = {}
        for svc, val in services_raw.items():
            if isinstance(val, (list, tuple)):
                if len(val) == 0:
                    logger.warning(
                        f"Service '{svc}' has an empty list; defaulting to (None, 1)"
                    )
                    sid, w = None, 1
                elif len(val) == 1:
                    sid, w = val[0], 1
                else:
                    sid, w = val[0], val[1]
            else:
                # scalar
                sid, w = val, 1

            normalized[svc] = (sid, w)
        return normalized

    ymllocation = os.path.join(downloaddir, "parseyml.yml")
    try:
        with open(ymllocation, "r") as stream:
            services_doc = yaml.safe_load(
                stream
            )  # equivalent to yaml.load(..., Loader=SafeLoader)

        if not isinstance(services_doc, dict) or "Services" not in services_doc:
            raise KeyError("'Services' key not found in YAML")

        ServicesRaw = services_doc["Services"]
        ServicesPairs = normalize_services(ServicesRaw)  # service -> (id, weight)
        Services = {svc: sid for svc, (sid, _) in ServicesPairs.items()}
        ServicesRaw = {
            svc: [sid, weight] for svc, (sid, weight) in ServicesPairs.items()
        }

    except FileNotFoundError:
        logger.error(f"YAML file not found: {ymllocation}", stack_info=True)
    except (yaml.YAMLError, KeyError, TypeError) as e:
        logger.error(f"Error reading or parsing YAML file: {e}", stack_info=True)

    ymllocation_hc = os.path.join(downloadhcdir, "parseyml.yml")

    with open(ymllocation_hc, "r") as stream:
        try:
            # Converts yaml document to python object
            services_hc = yaml.load(stream, Loader=SafeLoader)
            ServicesRaw_hc = services_hc["Services"]
            # Converts contents to usable dictionary
            Services_hc = {
                svc: (val[0] if isinstance(val, (list, tuple)) else val)
                for svc, val in ServicesRaw_hc.items()
            }
        except yaml.YAMLError as e:
            logger.error(f"Error reading YAML file: {e}", stack_info=True)
            raise e

    if args.years:
        config["hcstarty"] = args.years[0]
        config["hcendy"] = args.years[1]
    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016

    if args.yearsfc:
        years = [int(yr) for yr in args.yearsfc.split(",")]
        config["fcstarty"] = years[0]
        config["fcendy"] = years[0]
    else:
        config["fcstarty"] = 1993
        config["fcendy"] = 2016

    # hindcast info
    if centre == "eccc":
        # two models aka systems are live - call twice with each system number
        config["systemfc"] = Services["eccc_can"]
        config["systemhc"] = Services_hc["eccc_can"]
        compute_forecast(config, downloaddir, productshcdir, productsfcdir)

        ## repeat for second system
        config["systemfc"] = Services["eccc_gem5"]
        config["systemhc"] = Services_hc["eccc_gem5"]
        compute_forecast(config, downloaddir, productshcdir, productsfcdir)
    elif centre == "mme":
        config["systemfc"] = Services["mme"]
        ## Calculate fractional weights
        sum_weights = {}
        weights_sum = sum(n for _, n in ServicesRaw.values())
        for serv, val in ServicesRaw.items():
            version = val[0] if isinstance(val, (list, tuple)) else val
            weight = val[1] if isinstance(val, (list, tuple)) and len(val) > 1 else 1
            ServicesRaw[serv] = [version, weight / weights_sum]
        ## Run mme calc
        config["systemfc"] = Services["mme"]
        mme_products(ServicesRaw, config, productsfcdir)
    else:
        if centre not in Services.keys():
            logger.error(f"Unknown system for C3S: {centre}")
            raise ValueError(f"Unknown system for C3S: {centre}")
        config["systemfc"] = Services[centre]
        config["systemhc"] = Services_hc[centre]
        compute_forecast(config, downloaddir, productshcdir, productsfcdir)
    logger.info("Completed forecast products successfully")
