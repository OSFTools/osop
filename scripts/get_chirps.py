# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Python script to download CHIRPS data for a specific month, area and variable."""

import argparse
from datetime import datetime
import logging
from pathlib import Path

import requests

logger = logging.getLogger(__name__)


def get_obs(downloaddir, config):
    """Download CHIRPS for the requested period and area.

    Retrieves monthly averaged CHIRPS data. No server side subsetting is available,
    so the full dataset is downloaded and then subsetted to the requested area.

    Parameters
    ----------
    obs_fname : str
        Filename to save to, and to check not already downloaded.
    config : dict
        Dictionary containing necessary arguments.

    Returns
    -------
    str
        The filename where data is saved.
    """
    BASE_URL = "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/"

    starty = config["hcstarty"]
    endy = config["hcendy"]

    now = datetime.now()

    if starty < 1981:
        logger.warning(
            f"Data from before 1981 not available, setting start year to 1981"
        )
        starty = 1981
    if endy > now.year:
        logger.warning(
            f"Data from after {now.year} not available, setting end year to {now.year}"
        )
        endy = now.year

    nfail = 0
    for year in range(starty, endy + 1):
        for month in range(1, 13):
            logger.info(f"Downloading CHIRPS for {year}-{month:02d}")
            # Download the data for the specified year and month
            # Here we just simulate the download with a log message.
            obs_filename = f"chirps-v3.0.{year}.{month:02d}.tif"
            obs_fullpath = Path(downloaddir) / obs_filename

            if obs_fullpath.exists():
                message = f"File {obs_fullpath} already exists, skipping download."
                logger.warning(message)
                print(message)
                continue

            # now do the actual download (construct URL robustly to avoid double slashes)
            url = BASE_URL.rstrip("/") + "/" + obs_filename
            try:
                # use a short timeout and stream to avoid loading large responses into memory
                response = requests.get(url, timeout=30, stream=True)
                response.raise_for_status()  # Raise an error for bad responses
                with open(obs_fullpath, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                message = f"Successfully downloaded CHIRPS data for {year}-{month:02d} to {obs_fullpath}"
                logger.info(message)
                print(message)
            except requests.exceptions.RequestException as e:
                # log the URL and status if available to help debugging 403/forbidden
                status = (
                    getattr(e.response, "status_code", None)
                    if getattr(e, "response", None)
                    else None
                )
                message = f"Failed to download CHIRPS data for {year}-{month:02d}: {e} (url={url}, status={status})"
                logger.error(message)
                print(message)
                nfail += 1
    if nfail > 0:
        message = f"Finished downloading CHIRPS data with {nfail} failures."
        logger.error(message)
        raise RuntimeError(message)


def parse_args():
    """Set up argparse to get command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--month", required=True, help="start month for observations")
    parser.add_argument(
        "--leads", required=True, help="forecast range in months (comma separated)"
    )
    parser.add_argument(
        "--area",
        required=True,
        help="sub-area in degrees for retrieval (comma separated N,W,S,E)",
    )
    parser.add_argument("--downloaddir", required=True, help="location to download to")
    parser.add_argument("--logdir", required=True, help="location to store logfiles")
    parser.add_argument(
        "--pycptdir", required=True, help="location to store pycpt files"
    )
    parser.add_argument(
        "--pycpt", required=True, help="pycpt calibration: True or False"
    )
    parser.add_argument(
        "--years",
        required=False,
        help="Start and end years to retrieve data for (comma separated). Optional. Default is hindcast period 1993-2016.",
    )

    args = parser.parse_args()
    return args


def unpack_args_and_run(args):

    # unpack args and reformat if needed
    downloaddir = args.downloaddir
    pycptdir = args.pycptdir
    pycpt = args.pycpt

    # start logging - need to know logdir location before we can set it up
    logfile = (
        Path(args.logdir)
        / f"chirps_log_{datetime.today().strftime('%Y-%m-%d_%H:%M:%S')}.txt"
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

    month = int(args.month)
    leadtime_month = [int(l) - 1 for l in args.leads.split(",")]
    # for filename to keep consistent with hindcast filenames
    leads_str = "".join([str(mon) for mon in leadtime_month])
    area_bounds = [float(pt) for pt in args.area.split(",")]

    # add arguments to config dictionary used to pass parameters
    config = dict(
        start_month=month,
        leads_obs=leadtime_month,
        area=area_bounds,
        leads_str=leads_str,
    )

    logger.debug(config)

    if args.years:
        config["hcstarty"] = int(args.years.split(",")[0])
        config["hcendy"] = int(args.years.split(",")[1])

    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016

    get_obs(downloaddir, config)

    if pycpt == "True":
        raise NotImplementedError(
            "pycpt calibration not yet implemented for ERA5 retrievals"
        )


if __name__ == "__main__":
    """
    Called when this is run as a script.

    Gets the command line arguments using argparse and calls the main function to download ERA5.

    Returns
    -------
    None
    """
    # get command line args
    args = parse_args()
    unpack_args_and_run(args)
