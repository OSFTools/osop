"""
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

Python script to get NMME models not avialable from C3S
Using an interface as close to get_any_hc as possible

usage: get_any_hc-fc_iri.py.py [-h] --centre CENTRE --month MONTH --leads LEADS
                           --area AREA --downloaddir DOWNLOADDIR

options:
  -h, --help            show this help message and exit
  --model MODEL         model to download
  --month MONTH         start month for hindcasts/forecasts
  --leads LEADS         forecast range in months (comma separated)
  --area AREA           sub-area in degrees for retrieval (comma separated
                        N,W,S,E)
  --forecast            optional - default to false i.e. hindcast
  --downloaddir DOWNLOADDIR
                        location to download to
"""
import os
import sys
import datetime as dt
import cptdl
import argparse
import calendar

def round_month(month):
    """
    for month calculation, take remained when divided
    by 12 but for multiples of 12, return 12

    >>> print(round_month(1))
    1
    >>> print(round_month(12))
    12
    >>> print(round_month(13))
    1
    
    """
    month_out = month%12
    if month_out == 0:
         month_out = 12
     
    return month_out

def get_target(month, leadtime_month):
    """
    Given the month + leadtime_month variables, 
    calculate a string such as Dec-Feb for the target
    Using the same convention as for C3S forecasts 
    where month = 1 means the first month of the forecast
    
    >>> month = 11
    >>> leadtime_month =[2, 3, 4]
    >>> print(get_target(month, leadtime_month))
    Dec-Feb

    >>> month = 5
    >>> print(get_target(month, leadtime_month))
    Jun-Aug

    """
    
    month1 = round_month(month + leadtime_month[0] - 1)
    month2 = round_month(month + leadtime_month[-1] - 1)
    
    targ1 = calendar.month_name[month1][0:3]
    targ2 = calendar.month_name[month2][0:3]
    
    target = targ1 + '-' + targ2
    
    return target

def get_iri_nnme(model, month, leadtime_month, l_fc, area, area_str, downloaddir, year=None):

    VARIABLES = ["PRCP","T2M"]
    
    if l_fc:
        if year is None:
            raise ValueError('If downloading FC, need to specify year')
        start_time = dt.datetime(year,month,1)
    else:
        # Could be improved here - make dictionaries seperately
        start_time = dt.datetime(1,month,1)

    # string for the months to download
    leads_str = "".join([str(mon) for mon in leadtime_month])

	# kwargs common for all models (for forecasts)
    kwargs_fc = { 
        'fdate': start_time,
        'first_year': year, 
        'final_year': year, 
        # area is  N,W,S,E
        'predictor_extent': {
            'north': area[0],
            'west': area[1],
            'south': area[2],
            'east': area[3]
        }, 
        # for IRI, the month is the middle of the month of interest,
        # this is different to C3S where month 1 is the month the forecast
        # was initialised
        'lead_low': float(leadtime_month[0])-0.5,
        'lead_high': float(leadtime_month[-1])-0.5, 
        'target': get_target(month, leadtime_month),
        'filetype': 'cptv10.tsv',
        'ensemblemean': False
	}

	# kwargs for hindcasts - mostly the same
    kwargs_hc = kwargs_fc.copy()
	# key difference is that we want to get the 
	# data for the same period as from C3S
    # [perhaps alter later]
    kwargs_hc['first_year'] = 1993
    kwargs_hc['final_year'] = 2016
    #del kwargs_hc['fdate']

    
    # download one variable at a time as that is how IRI works
    for variable in VARIABLES:
        if l_fc:
            try:
                template = cptdl.forecasts[f'{model}.{variable}']
            except KeyError:
                print(f"Template not found for {model} - {variable}")
            
            # format of C3S download is {centre}_{system}_hc_monthly_mean_{month}_{leads_str}_{area_str}
            destination_file = f"{downloaddir}/{model}_fc_monthly_mean_{variable}_{year}{month}_{leads_str}_{area_str}.tsv" 
            print(f'downloading forecast: {destination_file}')

            # get as xarray dataarray
            da = cptdl.download(template, destination_file, use_dlauth=False, verbose=False, **kwargs_fc)  

            # want to include location downloaded from (base URL) in the global metadata
            base_url = cptdl.evaluate_url(template, **kwargs_hc)
            da.attrs['source'] = base_url

            # save as netCDF for further processing 
            da.to_netcdf(destination_file.replace('.tsv','.nc'))

            # delete tsv file
            os.remove(destination_file)
        else:
            try:
                template = cptdl.hindcasts[f'{model}.{variable}']
            except KeyError:
                print(f"Template not found for {model} - {variable}")

            # format of C3S download is {centre}_{system}_hc_monthly_mean_{month}_{leads_str}_{area_str}
            destination_file = f"{downloaddir}/{model}_hc_monthly_mean_{variable}_{month}_{kwargs_hc['first_year']}-{kwargs_hc['final_year']}_{leads_str}_{area_str}.tsv" 
            print(f'downloading hindcast: {destination_file}')
            print(cptdl.evaluate_url(template, **kwargs_hc)) 
            # get as xarray data array
            da = cptdl.download(template, destination_file, use_dlauth=False, verbose=False, **kwargs_hc)  

            # save as netCDF for further processing 
            # want to include location downloaded from (base URL) in the global metadata
            base_url = cptdl.evaluate_url(template, **kwargs_hc)
            da.attrs['source'] = base_url

            # ensemble member - don't store as string (crashes CF checker!)
            da["M"] = da.M.astype('float')
            da["M"] = da.M.astype('int')

            # also use metadata from cptdl
            for meta in ['institution', 'hindcast_limits', 'forecast_limits']:
                da.attrs[meta] = str(cptdl.catalog.metadata[model][meta])

            # change the start and end dates of the filename to reflect the actual contents
            y1 = str(da.S.data[0].astype('datetime64[Y]'))
            y2 = str(da.S.data[-1].astype('datetime64[Y]'))
            nc_file =  f"{downloaddir}/{model}_hc_monthly_mean_{variable}_{start_time:%d-%m}_{y1}-{y2}_{leads_str}_{area_str}.nc" 
            da.to_netcdf(nc_file)

            # delete tsv file
            os.remove(destination_file)

def parse_args():
    """
    set up argparse to get command line arguments
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, help="model to download")
    parser.add_argument("--month", required=True, help="start month for hindcasts")
    parser.add_argument(
        "--leads", required=True, help="forecast range in months (comma separated)"
    )
    parser.add_argument("--forecast", required=False, help="Get forecast? If false, get hindcast", action="store_true")
    # year only required if --forecast is given
    parser.add_argument('--year', required='--forecast' in sys.argv) 
    parser.add_argument(
        "--area",
        required=True,
        help="sub-area in degrees for retrieval (comma separated N,W,S,E)",
    )
    parser.add_argument("--downloaddir", required=True, help="location to download to")

    args = parser.parse_args()
    return args

def main():
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the main funciton to do the actual
    cptdl call to IRI
    """

    # get command line args
    args = parse_args()
    
    # unpack args and reformat if needed
    model = args.model
    downloaddir = args.downloaddir
    leadtime_month = [int(l) for l in args.leads.split(",")]
    area = [float(pt) for pt in args.area.split(",")]
    # create a string of the area to use in filenames
    area_str = args.area.replace(",", ":")
    month = int(args.month)

    if len(area) != 4:
        raise ValueError(f"Need 4 points for area: {area}")

    l_fc = args.forecast

    if l_fc:
        year = int(args.year)
    else:
        year = None

    get_iri_nnme(model, month, leadtime_month, l_fc, area, area_str, 
                downloaddir, year=year)
    
if __name__ == "__main__":
    main()

