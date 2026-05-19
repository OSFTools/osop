# (C) Crown Copyright, Met Office. All rights reserved.
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

import pycpt
import packaging
import numpy as np
from pathlib import Path
import xarray as xr
import eccodes

# my imports added below as needed
import pandas as pd

#For now starting with direct imports - can add function to pull through osop later 

ecmwf = '/data/scratch/eleanor.dean/seafoam/data/master/hindcast/pycpt/pycpt_ecmwf_51_1993-2016_monthly_mean_5_234_40:10:10:40_total_precipitation.pycpt.nc'
ukmo  = '/data/scratch/eleanor.dean/seafoam/data/master/hindcast/pycpt/pycpt_ukmo_604_1993-2016_monthly_mean_5_234_40:10:10:40_total_precipitation.pycpt.nc'
obs   = '/data/scratch/eleanor.dean/seafoam/data/master/hindcast/pycpt/era5_total_precipitation_1993-2016_monthly_5_123_40:10:10:40.pycpt.nc'

ecmwf_forecast = '/data/scratch/eleanor.dean/seafoam/data/master/forecast/pycpt/pycpt_ecmwf_51_2025-2025_monthly_mean_5_234_40:10:10:40_total_precipitation.pycpt.nc'
ukmo_forecast  = '/data/scratch/eleanor.dean/seafoam/data/master/forecast/pycpt/pycpt_ukmo_604_2025-2025_monthly_mean_5_234_40:10:10:40_total_precipitation.pycpt.nc'

#For now also forcing the case dir to just be a random one in sea foam - will correct later 

case_dir = Path('/data/scratch/eleanor.dean/seafoam/data/master/hindcast/pycpt')
print(case_dir)

#This domain will need a function to pull it through from OSOP
predictor_extent = {
    'west': 10,
    'east': 40,
    'south': 10,
    'north': 40,
}


#From here its forced set up of PYCPT settings - will need to discuss with NIck how I plan to put these in the master script. 
#Do we always want it turned on?

MOS = 'PCR'  # must be one of 'CCA', 'PCR'

cpt_args = {
    'transform_predictand': 'Gamma',
    'tailoring': None,
    'cca_modes': (1, 3),
    'x_eof_modes': (1, 8),
    'y_eof_modes': (1, 6),
    'validation': 'crossvalidation',
    'drymask_threshold': None,
    'skillmask_threshold': None,
    'crossvalidation_window': 5,
    'synchronous_predictors': True,
}

interactive = False




domain_dir = pycpt.setup(case_dir, predictor_extent)
pycpt.plot_domains(predictor_extent, predictor_extent)


#Loading in the local data - again will be function and generalised later to run in the OSOP script
ds_ecmwf = xr.open_dataset(ecmwf)
ds_ukmo  = xr.open_dataset(ukmo)
ds_obs   = xr.open_dataset(obs)

ds_ecmwf_fc = xr.open_dataset(ecmwf_forecast)
ds_ukmo_fc  = xr.open_dataset(ukmo_forecast)

# Extract variables (assumes one variable per file)
X_ecmwf = list(ds_ecmwf.data_vars.values())[0]
X_ukmo  = list(ds_ukmo.data_vars.values())[0]
Y       = list(ds_obs.data_vars.values())[0]

F_ecmwf = list(ds_ecmwf_fc.data_vars.values())[0]
F_ukmo  = list(ds_ukmo_fc.data_vars.values())[0]




# predictor_names provides the order, and hindcast_data / forecast_data must match it by index
predictor_names = ["ECMWF", "UKMO"]
predictand_name = "ERA5.PRCP"

hindcast_data = [X_ecmwf, X_ukmo]
forecast_data  = [F_ecmwf, F_ukmo]   # IMPORTANT: list, not dict, not None



hcsts, fcsts, skill, pxs, pys = pycpt.evaluate_models(
    hindcast_data,
    MOS,
    Y,
    forecast_data,
    cpt_args,
    domain_dir,
    predictor_names,
    interactive
)



skill_metrics = [
    "pearson",
    "spearman",
    "two_alternative_forced_choice",
    "roc_area_below_normal",
    "roc_area_above_normal",
    "root_mean_squared_error",
]

pycpt.plot_skill(predictor_names, skill, MOS, domain_dir, skill_metrics)
pycpt.plot_eof_modes(MOS, predictor_names, pxs, pys, domain_dir)
pycpt.plot_cca_modes(MOS, predictor_names, pxs, pys, domain_dir)
pycpt.plot_forecasts(cpt_args, predictand_name, fcsts, domain_dir, predictor_names, MOS)
ensemble = predictor_names

det_fcst, pr_fcst, pev_fcst, nextgen_skill = pycpt.construct_mme(fcsts, hcsts, Y, ensemble, predictor_names, cpt_args, domain_dir)

mme_skill_metrics = [
    "spearman",
    "two_alternative_forced_choice",
    "root_mean_squared_error",
    "generalized_roc",
    "rank_probability_skill_score",
]


pycpt.plot_mme_skill(predictor_names, nextgen_skill, MOS, domain_dir, mme_skill_metrics)
pycpt.plot_mme_forecasts(cpt_args, predictand_name, pr_fcst, MOS, domain_dir, det_fcst)

threshold = 0.5
isPercentile = True

exceedance_prob, fcst_scale, climo_scale, fcst_mu, climo_mu, Y2, ntrain, transformed_threshold = pycpt.construct_flex_fcst(MOS, cpt_args, det_fcst, threshold, isPercentile, Y, pev_fcst)

if 'station' in Y.coords:
    print(Y['Name'].to_pandas().to_string())

# location_selector = {'X': 42.75, 'Y': -75.88}  # example for gridded data
# location_selector = {'station': '11602'}  # example for station data
location_selector = None
forecast_colors = None

pycpt.plot_mme_flex_forecast_v2(predictand_name, exceedance_prob, transformed_threshold, fcst_scale, climo_scale, fcst_mu, climo_mu, Y, Y2, ntrain, MOS, domain_dir, location_selector=location_selector, color_bar=forecast_colors)

print("Run complete. Outputs saved to:", case_dir)