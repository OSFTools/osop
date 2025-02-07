"""
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

Script to download monthly mean seasonal data
from C3S for the common period of hindcasts 
1993-2016 or for any arbitrary period.

Usage:
python get_any_hindcast.py --centre <centre> --month <month> --leads <leads> 
                           --area <area> --downloaddir <downloaddir> --years <years>
    centre: modelling centre
    month: start month of forecast
    leads: forecast range in months (comma separated)
    area: sub-area in degrees for retrieval (comma separated N,W,S,E)
    downloaddir: location to download to
    years: Years to rerieve data for (comma separated). Optional. Default is hindcast period 1993-2016.

(C) Crown Copyright, Met Office. All rights reserved.
"""

import cdsapi
import argparse
import sys

sys.path.append('/home/h01/edyer/IASAS/osop/')
from lib.osop.constants import SYSTEMS


def do_cdsapi_call(
    centre, system, month, leadtime_month, area, area_str, variable, downloaddir, years="hc"
):
    """
    calls cdsapi for the requested period and area
    retrieves monthly mean seasonal data of t2m and total precip
    from the requested model

    Args:
        centre(str): modelling centre
        system(str): model version
        month(int): start month of forecast
        leadtime_month(list): list of lead times for FC (int)
        area(list): 4 element list of floats with area for download
        area_str(str): colon seperated area string
        downloaddir(str):directory to use for the downloads
        years(str): years to get data for (comma separated). Optional.
                    Default is hindcast period 1993-2016.


    TO DO:
    1. optionally, get other variables
    2. tests
    mocked:
    call with ecmwf and ecc
    different lengths of lead time

    should be failures:
    centre = 'bob'
    length not equal to 4

    """
    c = cdsapi.Client()

    leads_str = "".join([str(mon) for mon in leadtime_month])

    # set a default set of year that are the common period of hindcasts
    if years == "hc":
        years = [
            "1993",
            "1994",
            "1995",
            "1996",
            "1997",
            "1998",
            "1999",
            "2000",
            "2001",
            "2002",
            "2003",
            "2004",
            "2005",
            "2006",
            "2007",
            "2008",
            "2009",
            "2010",
            "2011",
            "2012",
            "2013",
            "2014",
            "2015",
            "2016",
        ]

    c.retrieve(
        "seasonal-monthly-single-levels",
        {
            "format": "grib",
            "originating_centre": f"{centre}",
            "system": system,
            "variable": [variable],
            "product_type": "monthly_mean",
            "year": years,
            "month": '{:02d}'.format(month),
            "leadtime_month": leadtime_month,
            "area": area,
        },
        f"{downloaddir}/{centre}_{system}_{years[0]}-{years[-1]}_monthly_mean_{month}_{leads_str}_{area_str}_{variable}.grib",
    )


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
        help="variable to download, 2m_temperature, total_precipitation")
    parser.add_argument("--downloaddir", required=True, help="location to download to")
    parser.add_argument(
        "--years",
        required=False,
        help="Years to rerieve data for (comma separated). Optional. Default is hindcast period 1993-2016.",
    )

    args = parser.parse_args()
    return args


def main():
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the main funciton to do the actual
    cdsapi call
    """

    # get command line args
    args = parse_args()

    # unpack args and reformat if needed
    centre = args.centre
    downloaddir = args.downloaddir
    leadtime_month = [int(l) for l in args.leads.split(",")]
    area = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    month = int(args.month)
    variable = str(args.variable)

    area = [float(pt) for pt in args.area.split(",")]
    if len(area) != 4:
        raise ValueError(f"Need 4 points for area: {area}")
    if args.years:
        years = [int(yr) for yr in args.years.split(",")]
    else:
        years = "hc"

    if centre == "eccc":
        # two models aka systems are live - call twice with each system number
        do_cdsapi_call(
            "eccc",
            SYSTEMS["eccc_can"],
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
            SYSTEMS["eccc_gem5"],
            month,
            leadtime_month,
            area,
            area_str,
            variable,
            downloaddir,
            years,
        )
    else:
        if centre not in SYSTEMS.keys():
            raise ValueError(f"Unknown system for C3S: {centre}")
        do_cdsapi_call(
            centre,
            SYSTEMS[centre],
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
