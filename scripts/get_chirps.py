# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Python script to download CHIRPS data for requested months and years."""

import argparse
from datetime import datetime
from datetime import timedelta
import logging
from pathlib import Path
import os

from dateutil.relativedelta import relativedelta
import requests
import rioxarray  # noqa: F401 - required to register rasterio engine with xarray
import xarray as xr
import numpy as np

from osop.pycpt_convert import process_grib_to_pycpt

logger = logging.getLogger(__name__)


def get_obs(obs_fname, config):
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

    if os.path.exists(obs_fname):
        logger.warning(f"File {obs_fname} already exists")
        return obs_fname

    # If hcstarty and hcendy exist in the config, use them to set the start and end years, else use 1993-2016
    try:
        starty = config["hcstarty"]
    except KeyError:
        starty = 1993
    starty0 = starty
    try:
        endy = config["hcendy"]
    except KeyError:
        endy = 2016
    endy0 = endy
    # starty = config["hcstarty"]
    # starty0 = starty
    # endy = config["hcendy"]
    # endy0 = endy

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
        datetime(iy, config["start_month"], 1) + relativedelta(months=ioff)
        for iy in range(starty, endy + 1)
        for ioff in config["leads_obs"]
    ]

    nfail = 0

    skipped = []
    succeeded = []
    subset_chirps_files = []

    print("Downloading CHIRPS data.")

    for dt in mon_dt:
        # Download the data for the specified year and month.
        year = dt.year
        month = dt.month
        logger.info(f"Downloading CHIRPS for {dt}")
        obs_filename = f"chirps-v3.0.{year}.{month:02d}.tif"
        obs_fullpath = Path(config["downloaddir"]) / obs_filename
        obs_nc_filename = f"chirps-v3.0.{year}.{month:02d}_f{config['area_str']}.nc"
        obs_fullpath_nc = Path(config["downloaddir"]) / obs_nc_filename

        if obs_fullpath.exists() or obs_fullpath_nc.exists():
            message = f"File {obs_fullpath_nc} already exists, skipping .tif download."
            print(message)
            logger.warning(message)
            skipped.append(obs_filename)

            sub_file = subset_chirps(
                obs_fullpath, config["area"], config["area_str"], config["ldelete"], year, month
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
                obs_fullpath, config["area"], config["area_str"], config["ldelete"], year, month
            )
            subset_chirps_files.append(sub_file)

        except requests.exceptions.RequestException as e:
            # Try to download the file again if it fails the first time
            message = f"Initial download failed for {year}-{month:02d}: {e}. Retrying with increased timeout..."
            logger.warning(message)
            print(message)
            try:
                # use a short timeout and stream to avoid loading large responses into memory
                response = requests.get(url, timeout=45, stream=True)
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
                    obs_fullpath, config["area"], config["area_str"], config["ldelete"], year, month
                )
                subset_chirps_files.append(sub_file)

            # If the retry also fails, log the error and increment the failure count.
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
  
    # Create one combined NetCDF file for all subsetted CHIRPS files
    combined_nc_file = obs_fname.replace(".grib", ".nc")
    
    if subset_chirps_files:
        combined_ds = xr.open_mfdataset(subset_chirps_files, combine="nested", concat_dim="start_date")
        # Convert negative tp values to NaN
        # combined_ds = combined_ds.where(combined_ds.tp >= 0)
        # Need to regrid to same grid used by C3S seasonal forecasts (1/1 degree)
        print("Regridding combined CHIRPS dataset.")
        combined_ds = regrid(combined_ds)
        combined_ds.to_netcdf(combined_nc_file)
        logger.info(f"Created combined subsetted CHIRPS file at {combined_nc_file}")
        print(f"Saved combined, subsetted and regridded CHIRPS file at: {combined_nc_file}")

def regrid(input_da):
    """Regrid the input DataArray to the standard 1/1 degree grid used by C3S seasonal forecasts."""
    from xarray_regrid import Grid

    # Get limits of dimensions from the input DataArray
    lon_min = input_da.longitude.min().round()
    lon_max = input_da.longitude.max().round()
    lat_min = input_da.latitude.min().round()
    lat_max = input_da.latitude.max().round()
    longitude = np.arange(lon_min, lon_max + 1, 1)
    latitude = np.arange(lat_min, lat_max + 1, 1)
    target_grid = Grid(
    north=lat_max,
    south=lat_min,
    east=lon_max,
    west=lon_min,
    resolution_lat=1,
    resolution_lon=1,
    )
    target_grid = target_grid.create_regridding_dataset()
    da_regridded = input_da.regrid.conservative(target_grid, time_dim="start_date")
   
    return da_regridded

def subset_chirps(tif_file, area_bounds, area_str, ldelete, year, month):
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
    year : int
        Year of the CHIRPS data.
    month : int
        Month of the CHIRPS data

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

    # # Open the GeoTIFF file using xarray with rasterio engine
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

    # Process subset to align attrs and units with ERA5
    subset = update_attrs(subset, year, month)
    
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

def update_attrs(chirps_file, year, month):
    """Update the attributes of CHIRPS data to align with ERA5 conventions for onward compatibility.

    Parameters
    ----------
    chirps_file : xarray.DataArray
        The CHIRPS data.
    year : int
        Year of the CHIRPS data.
    month : int
        Month of the CHIRPS data.

    Returns
    -------
    xarray.DataArray
        The subsetted CHIRPS data with updated attributes.
    """

    # Dictionary of attributes to update
    attr_dict = {
        "Conventions": "CF-1.6",
        "title": "CHIRPS Version 3.0",
        "history": "created by Climate Hazards Center",
        "version": "Version 3.0",
        "creator_name": "Pete Peterson",
        "creator_email": "pete@geog.ucsb.edu",
        "institution": "Climate Hazards Center. University of California at Santa Barbara",
        "documentation": "http://pubs.usgs.gov/ds/832/",
        "reference": "Funk, C.C., Peterson, P.J., Landsfeld, M.F., Pedreros, D.H., Verdin, J.P., Rowland, J.D., Romero, B.E., Husak, G.J., Michaelsen, J.C., and Verdin, A.P., 2014, A quasi-global precipitation time series for drought monitoring: U.S. Geological Survey Data Series 832, 4 p., http://dx.doi.org/110.3133/ds832. ",
        "comments":  "time variable denotes the first day of the given month. ",
        "acknowledgements": "The Climate Hazards Center InfraRed Precipitation with Stations development process was carreied out through U.S. Geological Survey (USGS) cooperative agreement #G09AC000001 'Monitoring and Forecasting Climate, Water and Land Use for Food Production in the Developing World' with funding from: U.S. Agency for International Development Office of Food for Peace, award #AID-FFP-P-10-00002 for 'Famine Early Warning Systems Network Support,' the National Aeronautics and Space Administration Applied Sciences Program, Decisions award #NN10AN26I for 'A Land Data Assimilation System for Famine Early Warning,' SERVIR award #NNH12AU22I for 'A Long Time-Series Indicator of Agricultural Drought for the Greater Horn of Africa,' The National Oceanic and Atmospheric Administration award NA110AR431051 for 'A Global Standardized Precipitation Index supporting the US Drought Portal and the Famine Early Warning System Network,' and the USGS Land Change Science Program.",
        "ftp_url": "ftp://chg-ftpout.geog.ucsb.edu/pub/chg/products/CHIRPS-latest/",
        "website": "http://chg.geog.ucsb.edu/data/chirps/index.html",
        "faq": "http://chg-wiki.geog.ucsb.edu/wiki/CHIRPS_FAQ",
    }

    tp_attr_dict = {
        "units": "m",
        "standard_name": "convective precipitation rate",
        "long_name": "Climate Hazards center InfraRed Precipitation with Stations",
        "time_step": "month",
        "geostatial_lat_min": -60,
        "geostatial_lat_max": 60,
        "geostatial_lon_min": -180,
        "geostatial_lon_max": 180,
    }
    
    # Get values from chirps_file before creating a new xr dataset
    chirps_tp_values = chirps_file.values
    new_shape = chirps_tp_values.shape[1:3]
    chirps_tp_values = chirps_tp_values.reshape(new_shape)
    chirps_lon_values = chirps_file.x.values
    chirps_lat_values = chirps_file.y.values

    chirps_coords = {"longitude": chirps_lon_values, "latitude": chirps_lat_values}
    chirps_var_name = "tp"
    chirps_ds = xr.Dataset(
        {chirps_var_name: (["latitude", "longitude"], chirps_tp_values)},
        coords=chirps_coords
    )

    # Update the attributes
    chirps_ds.attrs.update(attr_dict)

    chirps_ds = chirps_ds / 1000.0
    chirps_ds.tp.attrs = {}
    chirps_ds.tp.attrs.update(tp_attr_dict)

    # Update coordinate names and attributes
    chirps_ds.longitude.attrs = {'units': 'degrees_east', 'standard_name': 'longitude', 'long_name': 'longitude', 'axis': 'X'}
    chirps_ds.latitude.attrs = {'units': 'degrees_north', 'standard_name': 'latitude', 'long_name': 'latitude', 'axis': 'Y'}

    # Add valid_time with format YYYY-MM-01T06:00:00Z
    valid_time_val = f"{year}-{month:02d}-01T06:00:00Z"
    valid_time_val = datetime.strptime(valid_time_val, "%Y-%m-%dT%H:%M:%SZ")
    chirps_ds = chirps_ds.assign_coords(valid_time=valid_time_val)

    # Create start_date values in line with era5 data, these are valid_time - 12 hours
    start_date_val = valid_time_val - timedelta(hours=12)
    chirps_ds = chirps_ds.assign_coords(start_date=start_date_val)
    
    return chirps_ds

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

    parser.add_argument(
        "--variable",
        required=False,
        help="variable name expected from config, although only tp is available",
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
        / f"chirps_log_{datetime.today().strftime('%Y.%m.%d_%H.%M.%S')}.txt"
    )

    loglev = logging.INFO  # can be an argument later if needed
    logging.basicConfig(
        level=loglev,
        filename=logfile,
        encoding="utf-8",
        filemode="w",
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y.%m.%d %H.%M",
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

    # unpack args and reformat if needed
    downloaddir = args.downloaddir
    pycptdir = args.pycptdir
    pycpt = args.pycpt

    # start logging - need to know logdir location before we can set it up
    logfile = os.path.join(
        args.logdir,
        f"chirps_log_{args.variable}_{args.month}_{datetime.today().strftime('%Y.%m.%d_%H.%M.%S')}.txt",
    )
    loglev = logging.INFO  # can be an argument later if needed
    logging.basicConfig(
        level=loglev,
        filename=logfile,
        encoding="utf-8",
        filemode="w",
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y.%m.%d %H.%M",
    )

    month = int(args.month)
    leadtime_month = [int(l) - 1 for l in args.leads.split(",")]
    # for filename to keep consistent with hindcast filenames
    leads_str = "".join([str(mon) for mon in leadtime_month])
    area_bounds = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ".")
    var = args.variable

    # add arguments to config dictionary used to pass parameters
    config = dict(
        start_month=month,
        leads_obs=leadtime_month,
        area=area_bounds,
        area_str=area_str,
        leads_str=leads_str,
        var=var,
        downloaddir=downloaddir,
        ldelete=True,
    )

    logger.debug(config)

    if args.years:
        config["hcstarty"] = int(args.years[0])
        config["hcendy"] = int(args.years[1])

    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016

    ## obs info
    obs_fname = "{fpath}/chirps_{var}_{hcstarty}-{hcendy}_monthly_{start_month}_{leads_str}_{area_str}.grib".format(
        fpath=downloaddir, **config
    )

    logger.info(f"Downloading obs filename: {obs_fname}")
    get_obs(obs_fname, config)

    if pycpt == "True":
        process_grib_to_pycpt(
            config,
            downloaddir,
            pycptdir,
            "obs",
            steps_to_sum=len(config["leads_obs"]),
            lead_months=1,
        )




    # unpack_args_and_run(args)
