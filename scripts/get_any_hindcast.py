# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Script to download monthly mean seasonal data.

Download from C3S for the common period of hindcasts 1993-2016 or for any arbitrary period.

Usage:
    python get_any_hindcast.py --centre <centre> --month <month> --leads <leads> \
        --area <area> --downloaddir <downloaddir> --years <years>
    centre: modelling centre
    month: start month of forecast
    leads: forecast range in months (comma separated)
    area: sub-area in degrees for retrieval (comma separated N,W,S,E)
    downloaddir: location to download to
    years: Years to retrieve data for (comma separated). Optional. Default is hindcast period 1993-2016.
"""

import argparse
from datetime import datetime
import logging
import os

import cdsapi
import yaml
from yaml.loader import SafeLoader

# Ensure the top level directory has been added to PYTHONPATH

logger = logging.getLogger(__name__)


def do_cdsapi_call(
    centre,
    system,
    month,
    leadtime_month,
    area,
    area_str,
    variable,
    downloaddir,
    years="hc",
):
    """Call cdsapi for the requested period and area.

    Retrieve monthly mean seasonal data of t2m and total precip from the requested model.

    Args:
        centre(str): modelling centre
        system(str): model version
        month(int): start month of forecast
        leadtime_month(list): list of lead times for FC (int)
        area(list): 4 element list of floats with area for download
        area_str(str): colon seperated area string
        variable(str): variable for download
        downloaddir(str):directory to use for the downloads
        years(str): years to get data for (comma separated). Optional.
                    Default is hindcast period 1993-2016.

    """
    leads_str = "".join([str(mon) for mon in leadtime_month])

    # set a default set of years that are the common period of hindcasts
    if years == "hc":
        years = [str(y) for y in range(1993, 2017)]

    fname = f"{downloaddir}/{centre}_{system}_{years[0]}-{years[-1]}_monthly_mean_{month}_{leads_str}_{area_str}_{variable}.grib"
    if os.path.exists(fname):
        logger.warning(f"File {fname} already exists")

    else:
        c = cdsapi.Client()
        c.retrieve(
            "seasonal-monthly-single-levels",
            {
                "format": "grib",
                "originating_centre": f"{centre}",
                "system": system,
                "variable": [variable],
                "product_type": "monthly_mean",
                "year": years,
                "month": "{:02d}".format(month),
                "leadtime_month": leadtime_month,
                "area": area,
            },
            f"{downloaddir}/{centre}_{system}_{years[0]}-{years[-1]}_monthly_mean_{month}_{leads_str}_{area_str}_{variable}.grib",
        )
        logger.info(f"Downloaded {fname}")


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
    parser.add_argument("--downloaddir", required=True, help="location to download to")

    parser.add_argument(
        "--logdir",
        required=True,
        help="Path to log directory for logging output.",
    )

    parser.add_argument(
        "--years",
        required=False,
        help="Years to rerieve data for (comma separated). Optional. Default is hindcast period 1993-2016.",
    )

    args = parser.parse_args()
    return args


def main():
    """Execute the script to download seasonal forecast data.

    Get the command line arguments using argparse and call the cdsapi
    to download the requested data.
    """
    # get command line args
    args = parse_args()

    # start logging - need to know logdir location before we can set it up
    logfile = os.path.join(
        args.logdir,
        f"download_log_{args.variable}_{args.centre}_{args.month}_{datetime.today().strftime('%Y-%m-%d_%H:%M:%S')}.txt",
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
    leadtime_month = [int(l) for l in args.leads.split(",")]
    area = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    month = int(args.month)
    variable = str(args.variable)

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

    area = [float(pt) for pt in args.area.split(",")]
    if len(area) != 4:
        logger.error(f"Need 4 points for area: {area}", stack_info=True)
        raise ValueError(f"Need 4 points for area: {area}")
    if args.years:
        years = [int(yr) for yr in args.years.split(",")]
    else:
        years = "hc"

    if centre == "eccc":
        # two models aka systems are live - call twice with each system number
        do_cdsapi_call(
            "eccc",
            Services["eccc_can"],
            month,
            leadtime_month,
            area,
            area_str,
            variable,
            downloaddir,
            years,
        )
        do_cdsapi_call(
            "eccc",
            Services["eccc_gem5"],
            month,
            leadtime_month,
            area,
            area_str,
            variable,
            downloaddir,
            years,
        )
    else:
        if centre not in Services.keys():
            logger.error(f"Unknown system for C3S: {centre}", stack_info=True)
            raise ValueError(f"Unknown system for C3S: {centre}")
        do_cdsapi_call(
            centre,
            Services[centre],
            month,
            leadtime_month,
            area,
            area_str,
            variable,
            downloaddir,
            years,
        )


if __name__ == "__main__":
    main()
