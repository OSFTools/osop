"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.

Script to download monthly mean seasonal data 
from C3S for the common period of hindcasts 
1993-2016 or for any arbitrary period.

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

import cdsapi
import os
import yaml
from yaml.loader import SafeLoader

# Ensure the top level directory has been added to PYTHONPATH


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
        variable(str): variable for download
        downloaddir(str):directory to use for the downloads
        years(str): years to get data for (comma separated). Optional.
                    Default is hindcast period 1993-2016.

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
    fname = f"{downloaddir}/{centre}_{system}_{years[0]}-{years[-1]}_monthly_mean_{month}_{leads_str}_{area_str}_{variable}.grib"
    if os.path.exists(fname):
        print(f"File {fname} already exists")
    else:
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


def main():
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the main funciton to do the actual
    cdsapi call
    """
    # get remaning arguments from yml file
    ymllocation = os.path.join("variables.yml")

    with open(ymllocation, "r") as stream:
        try:
            print('yml found')
            # Converts yaml document to python object
            config_test = yaml.load(stream, Loader=SafeLoader)
            # Converts contents to useable dictionary
            Services = config_test["Services"]
            Month_test = config_test["month"]
            leads_test = config_test["leads"]
            area_test = config_test["area"]
            varaible_test = config_test["variable"]
            downloaddir_test = config_test["downloaddir"]
            years = config_test["years"]
            centre = config_test["centre"]
            print('yml success')
        except yaml.YAMLError as e:
            print(e)
    
    # unpack args and reformat if needed
    centre = centre
    downloaddir_test = os.path.expandvars(downloaddir_test)
    downloaddir = downloaddir_test
    leadtime_month = [int(l) for l in leads_test.split(",")]
    area = [float(pt) for pt in area_test.split(",")]
    area_str = area_test.replace(",", ":")
    month = int(Month_test)
    variable = str(varaible_test)

    area = [float(pt) for pt in area_test.split(",")]
    if len(area) != 4:
        raise ValueError(f"Need 4 points for area: {area}")
    if years:
        years = [int(yr) for yr in years.split(",")]
    else:
        years = "hc"
    
    for centre in config_test["centre"]:
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
