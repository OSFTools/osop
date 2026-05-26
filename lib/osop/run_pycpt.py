# (C) Crown Copyright, Met Office. All rights reserved.
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Runs the pycpt workflow when PYCPT=TRUE in the .sh scripts.

The workflow is currently static until review.
It takes the pycpt files formatted from pycpt_convert as inputs and outputs the files from the designated workflow.

"""

from pathlib import Path

import eccodes
import mpl
import numpy as np
import packaging

# my imports added below as needed
import pandas as pd
import pycpt
import xarray as xr

mpl.use("Agg")


def process_pycpt(
    predict_config,
    pycptdir,
    hindcast_pycptdir,
):
    """Open and create the variables/datasets for pycpt.

    predict_config(dict): Config for file names.
    pycptdir(str): Forecast directory.
    hindcast_pycptdir(str): Hindcast directory.

    Returns: None.
    """
    print(predict_config)

    # Call the hindcast pycpt file
    Hcast_bname = "{pycptver}_{origin}_{systemhc}_1993-2016_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.pycpt".format(
        **predict_config
    )
    Hcast_fname = f"{hindcast_pycptdir}/{Hcast_bname}.nc"

    # Call the forecast pycpt file
    Fcast_bname = "{pycptver}_{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.pycpt".format(
        **predict_config
    )
    Fcast_fname = f"{pycptdir}/{Fcast_bname}.nc"

    # Call the predictor file
    Ocast_bname = "era5_{hc_var}_1993-2016_monthly_{start_month}_{obs_str}_{obs_area}.pycpt".format(
        **predict_config
    )
    Ocast_fname = f"{hindcast_pycptdir}/{Ocast_bname}.nc"

    # Set up output dir
    output_case_dir = Path(hindcast_pycptdir)

    # Set up bounds
    north, west, south, east = map(float, predict_config["obs_area"].split(":"))

    predictor_extent = {
        "west": west,
        "east": east,
        "south": south,
        "north": north,
    }

    origin = predict_config["origin"].upper()
    predictor_names = [origin]

    var = predict_config["var"]
    if var == "tprate":
        obs_var = "PRCP"
    elif var == "t2m":
        obs_var = "T2M"
    else:
        raise ValueError(f"Unsupported var: {var}")

    predictand_name = f"ERA5.{obs_var}"

    calibrate(
        Hcast_fname,
        Fcast_fname,
        Ocast_fname,
        output_case_dir,
        predictor_extent,
        predictor_names,
        predictand_name,
    )
    return


def calibrate(
    Hfname, Ffname, Ofname, case_dir, predictor_extent, predictor_name, predictand_name
):
    """Run Pycpt on loaded datafiles.

    Hfname(str): Hindcast location and filename.
    Ffname(str): Forecast location and filename.
    Ofname(str): Obs location and filename.
    ouput_case_dir(str): Output dir.
    predictor_extent(str): Predictor area bounds.
    predictor_name(arr): Name of predictor.
    predictand_name(arr): Name of predictand.

    Returns: None.
    """
    # Loading in the local data - again will be function and generalised later to run in the OSOP script
    ds_hindcast = xr.open_dataset(Hfname)
    ds_forecast = xr.open_dataset(Ffname)
    ds_obs = xr.open_dataset(Ofname)

    # Extract variables (assumes one variable per file)
    X_hindcast = list(ds_hindcast.data_vars.values())[0]
    Y = list(ds_obs.data_vars.values())[0]
    F_forecast = list(ds_forecast.data_vars.values())[0]

    MOS = "PCR"  # must be one of 'CCA', 'PCR'

    cpt_args = {
        "transform_predictand": None,
        "tailoring": None,
        "cca_modes": (1, 3),
        "x_eof_modes": (1, 8),
        "y_eof_modes": (1, 6),
        "validation": "crossvalidation",
        "drymask_threshold": None,
        "skillmask_threshold": None,
        "crossvalidation_window": 5,
        "synchronous_predictors": True,
    }

    interactive = False
    domain_dir = pycpt.setup(case_dir, predictor_extent)
    # pycpt.plot_domains(predictor_extent, predictor_extent)

    # predictor_names provides the order, and hindcast_data / forecast_data must match it by index
    predictor_names = predictor_name
    predictand_name = predictand_name

    hindcast_data = [X_hindcast]
    forecast_data = [F_forecast]  # IMPORTANT: list, not dict, not None

    hcsts, fcsts, skill, pxs, pys = pycpt.evaluate_models(
        hindcast_data,
        MOS,
        Y,
        forecast_data,
        cpt_args,
        domain_dir,
        predictor_names,
        interactive,
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
    pycpt.plot_forecasts(
        cpt_args, predictand_name, fcsts, domain_dir, predictor_names, MOS
    )

    return
