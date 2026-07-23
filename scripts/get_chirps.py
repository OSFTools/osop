# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Python script to download CHIRPS data for requested months and years."""

import argparse
from datetime import datetime
import logging
from pathlib import Path

from dateutil.relativedelta import relativedelta
import requests
import rioxarray  # noqa: F401 - required to register rasterio engine with xarray
import xarray as xr

logger = logging.getLogger(__name__)


def get_obs(config):
    """Download CHIRPS for the requested period and area.

    Retrieves monthly averaged CHIRPS data. No server side subsetting is available,
    so the full dataset is downloaded (no subsetting yet).

    Parameters
    ----------
    downloaddir : str | pathlib.Path
        Directory to save downloaded global GeoTIFFs.
    config : dict
        Dictionary containing necessary arguments.

    Returns
    -------
    None
        Downloads files into downloaddir; raises RuntimeError if any downloads fail.
    """
    BASE_URL = "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/"

    starty = config["hcstarty"]
    starty0 = starty
    endy = config["hcendy"]
    endy0 = endy

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

    if starty != starty0 or endy != endy0:
        print(
            f"Adjusted requested years from {starty0}-{endy0} to {starty}-{endy} based on available data."
        )

    # set up a list of datetimes for the requested months and years

    mon_dt = [
        datetime(iy, config["month"], 1) + relativedelta(months=ioff)
        for iy in range(starty, endy + 1)
        for ioff in config["leadtime_month"]
    ]
    nfail = 0

    skipped = []
    succeeded = []
    subset_chirps_files = []
    for dt in mon_dt:
        # Download the data for the specified year and month.
        year = dt.year
        month = dt.month
        logger.info(f"Downloading CHIRPS for {dt}")
        obs_filename = f"chirps-v3.0.{year}.{month:02d}.tif"
        obs_fullpath = Path(config["downloaddir"]) / obs_filename

        if obs_fullpath.exists():
            message = f"File {obs_fullpath} already exists, skipping download."
            logger.warning(message)
            skipped.append(obs_filename)

            sub_file = subset_chirps(
                obs_fullpath, config["area"], config["area_str"], config["ldelete"]
            )
            subset_chirps_files.append(sub_file)
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
            succeeded.append(obs_filename)
            # subset the downloaded file to the requested area
            sub_file = subset_chirps(
                obs_fullpath, config["area"], config["area_str"], config["ldelete"]
            )
            subset_chirps_files.append(sub_file)

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

    print(
        f"Finished downloading CHIRPS data. {len(succeeded)} files downloaded, {len(skipped)} files skipped."
    )
    print(f"Downloaded files: {', '.join(succeeded)}")
    print(f"Skipped files: {', '.join(skipped)}")

    print(f"Subsetted files: {', '.join(subset_chirps_files)}")


def subset_chirps(tif_file, area_bounds, area_str, ldelete):
    """Subset a downloaded CHIRPS GeoTIFF to the specified area and write a NetCDF.

    Parameters
    ----------
    tif_file : pathlib.Path
        Path to the downloaded CHIRPS GeoTIFF.
    area_bounds : list[float]
        Bounds of the area to subset, in degrees (N, W, S, E).
    area_str : str
        Bounds string used in the output filename (e.g. "10.-20.-10.20").
    ldelete : bool
        If True, delete the source GeoTIFF after successfully writing the NetCDF.

    Returns
    -------
    str
        Path to the subset NetCDF file.
    """
    max_lat, min_lon, min_lat, max_lon = area_bounds

    # Derive output file path early so we can skip work if it already exists.
    nc_output_file = tif_file.parent / f"{tif_file.stem}_f{area_str}.nc"
    if nc_output_file.exists() and nc_output_file.stat().st_size > 0:
        logger.info(
            f"Subsetted file {nc_output_file} already exists, skipping subsetting."
        )
        return str(nc_output_file)

    # Open the GeoTIFF file using xarray with rasterio engine
    ds = xr.open_dataarray(str(tif_file), engine="rasterio")

    # Check longitude bounds and raise ValueError if they are outside the CHIRPS extent.
    x_min = ds.x.min()
    x_max = ds.x.max()
    if not (x_min <= min_lon <= x_max and x_min <= max_lon <= x_max):
        raise ValueError("Longitude bounds are out of range of the CHIRPS data.")
    # Subset the data to the specified area with a 1 point buffer to avoid data loss at the edges.
    # The buffer is applied by extending the slice by 1 index in each direction.
    dx = abs(ds.x[1] - ds.x[0])
    dy = abs(ds.y[1] - ds.y[0])
    subset = ds.sel(
        y=slice(max_lat + dy, min_lat - dy), x=slice(min_lon - dx, max_lon + dx)
    )

    # Save the subsetted data back to the new file
    subset.to_netcdf(nc_output_file)
    logger.info(f"Successfully opened and subsetted {tif_file} to {nc_output_file}")

    # delete the original tif file if needed, or keep it for reference
    if ldelete:
        try:
            tif_file.unlink()
            logger.info(f"Deleted original file {tif_file} after subsetting.")
        except OSError as e:
            logger.error(f"Failed to delete original file {tif_file}: {e}")

    return str(nc_output_file)


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

    parser.add_argument(
        "--lkeep",
        action="store_true",
        required=False,
        help="Keep original tif files after subsetting",
    )
    args = parser.parse_args()
    return args


def unpack_args_and_run(args):
    """Unpack command line arguments and call the main function to download CHIRPS."""
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
    area_bounds = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ".")
    if args.lkeep:
        ldelete = False
    else:
        ldelete = True

    # add arguments to config dictionary used to pass parameters
    config = dict(
        month=month,
        area=area_bounds,
        area_str=area_str,
        leadtime_month=leadtime_month,
        downloaddir=downloaddir,
        ldelete=ldelete,
    )

    logger.debug(config)

    if args.years:
        config["hcstarty"] = int(args.years.split(",")[0])
        config["hcendy"] = int(args.years.split(",")[1])

    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016
    get_obs(config)

    if pycpt == "True":
        raise NotImplementedError(
            "pycpt calibration not yet implemented for ERA5 retrievals"
        )


if __name__ == "__main__":
    """
    Called when this is run as a script.

    Gets the command line arguments using argparse and calls the main function to download CHIRPS.

    Returns
    -------
    None
    """
    # get command line args
    args = parse_args()
    unpack_args_and_run(args)
