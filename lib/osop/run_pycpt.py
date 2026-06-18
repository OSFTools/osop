# (C) Crown Copyright, Met Office. All rights reserved.
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Runs the pycpt workflow when PYCPT=TRUE in the .sh scripts.

It takes the pycpt files formatted from pycpt_convert as inputs and outputs the files from the designated workflow.

"""

import logging
from pathlib import Path
import warnings

logger = logging.getLogger(__name__)

import eccodes
import matplotlib as mpl
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
    ServicesRaw,
):
    """Open and create the variables/datasets for pycpt.

    predict_config(dict): Config for file names.
    pycptdir(str): Forecast directory.
    hindcast_pycptdir(str): Hindcast directory.

    Returns: None.
    """
    # Set up output dir
    output_case_dir = Path(hindcast_pycptdir)

    # Set up bounds (obs_area stored "N:W:S:E")
    north, west, south, east = map(float, predict_config["obs_area_str"].split(":"))
    predictand_extent = {
        "west": west,
        "east": east,
        "south": south,
        "north": north,
    }

    # Set up bounds for predictor
    north_M, west_M, south_M, east_M = map(
        float, predict_config["gcm_area_str"].split(":")
    )
    predictor_extent = {
        "west": west_M,
        "east": east_M,
        "south": south_M,
        "north": north_M,
    }

    var = predict_config["var"]
    if var == "tprate":
        obs_var = "PRCP"
    elif var == "t2m":
        obs_var = "T2M"
    else:
        raise ValueError(f"Unsupported var: {var}")

    predictand_name = f"ERA5.{obs_var}"

    # Call the predictor (Obs)
    Ocast_bname = "era5_{hc_var}_1993-2016_monthly_{start_month}_{obs_str}_{obs_area_str}.pycpt".format(
        **predict_config
    )
    Ocast_fname = f"{hindcast_pycptdir}/{Ocast_bname}.nc"
    if not Path(Ocast_fname).is_file():
        raise FileNotFoundError(f"Obs pycpt file not found: {Ocast_fname}")

    if not isinstance(ServicesRaw, dict):
        raise TypeError(
            "ServicesRaw must be a dict of service, ID and weight [id, weight]"
        )

    services_values = {
        svc: (val[0] if isinstance(val, (list, tuple)) else val)
        for svc, val in ServicesRaw.items()
    }
    services_weights = {
        svc: (val[1] if isinstance(val, (list, tuple)) and len(val) > 1 else 1)
        for svc, val in ServicesRaw.items()
    }

    # Skip any weight==0
    active_services = [
        svc for svc in services_values.keys() if services_weights.get(svc, 1) > 0
    ]

    hindcast_files = []
    forecast_files = []
    predictor_names = []

    for svc in active_services:
        # Skip MME.
        if svc == "mme":
            continue

        cfg = dict(predict_config)  # shallow copy is fine for string/int fields
        cfg["origin"] = svc

        # Use service id for file naming for both hc and fc unless already set and matching svc
        sid = services_values[svc]
        cfg["systemhc"] = sid
        cfg["systemfc"] = sid

        # Hindcast filename
        Hcast_bname = "{pycptver}_{origin}_{systemhc}_1993-2016_monthly_mean_{start_month}_{leads_str}_{gcm_area_str}_{hc_var}.pycpt".format(
            **cfg
        )
        Hcast_fname = f"{hindcast_pycptdir}/{Hcast_bname}.nc"

        # Forecast filename
        Fcast_bname = "{pycptver}_{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{gcm_area_str}_{hc_var}.pycpt".format(
            **cfg
        )
        Fcast_fname = f"{pycptdir}/{Fcast_bname}.nc"

        # Skip if either file missing
        if not Path(Hcast_fname).is_file():
            warnings.warn(
                f"Missing hindcast pycpt file for service '{svc}': {Hcast_fname} - skipping"
            )
            continue
        if not Path(Fcast_fname).is_file():
            warnings.warn(
                f"Missing forecast pycpt file for service '{svc}': {Fcast_fname} - skipping"
            )
            continue

        hindcast_files.append(Hcast_fname)
        forecast_files.append(Fcast_fname)
        predictor_names.append(svc.upper())

    if len(hindcast_files) == 0:
        raise ValueError(
            "No valid service pairs found (hindcast+forecast). "
            "Check your pycpt files exist in both directories, or update master.forecast.sh to generate them."
        )

    calibrate(
        hindcast_files,
        forecast_files,
        Ocast_fname,
        output_case_dir,
        predictor_extent,
        predictand_extent,
        predictor_names,
        predictand_name,
    )
    return


def calibrate(
    Hfnames,
    Ffnames,
    Ofname,
    case_dir,
    predictor_extent,
    predictand_extent,
    predictor_names,
    predictand_name,
):
    """Run Pycpt on loaded datafiles.

    Hfnames(list[str]): Hindcast locations and filenames (one per predictor).
    Ffnames(list[str]): Forecast locations and filenames (one per predictor).
    Ofname(str): Obs location and filename.
    ouput_case_dir(str): Output dir.
    predictor_extent(dict): Predictor area bounds.
    predictand_extent(dict): Predictand area bounds.
    predictor_names(list[str]): Names of predictors (must align with Hfnames/Ffnames order).
    predictand_name(str): Name of predictand.

    Returns: None.
    """
    # Load obs
    ds_obs = xr.open_dataset(Ofname)
    Y = list(ds_obs.data_vars.values())[0]

    hindcast_data = []
    forecast_data = []
    kept_names = []

    # Load each predictor pair
    for Hfname, Ffname, pname in zip(Hfnames, Ffnames, predictor_names):
        try:
            ds_hindcast = xr.open_dataset(Hfname)
            ds_forecast = xr.open_dataset(Ffname)

            # Extract variables (assumes one variable per file)
            X_hindcast = list(ds_hindcast.data_vars.values())[0]
            F_forecast = list(ds_forecast.data_vars.values())[0]

            hindcast_data.append(X_hindcast)
            forecast_data.append(F_forecast)
            kept_names.append(pname)
        except Exception as e:
            warnings.warn(
                f"Failed to load predictor '{pname}' (H='{Hfname}', F='{Ffname}'): {e} - skipping"
            )
            continue

    predictor_names = kept_names

    if len(hindcast_data) == 0:
        raise ValueError("No predictors could be loaded successfully.")

    MOS = "CCA"  # must be one of 'CCA', 'PCR'

    cpt_args = {
        "transform_predictand": None,
        "tailoring": None,
        "cca_modes": (1, 10),
        "x_eof_modes": (1, 10),
        "y_eof_modes": (1, 10),
        "validation": "crossvalidation",
        "drymask_threshold": None,
        "skillmask_threshold": None,
        "crossvalidation_window": 5,
        "synchronous_predictors": True,
    }

    interactive = False
    domain_dir = pycpt.setup(case_dir, predictor_extent)

    hcsts, fcsts, skill, pxs, pys = pycpt.evaluate_models(
        hindcast_data,
        MOS,
        Y,
        forecast_data,
        cpt_args,
        domain_dir,
        predictor_names,  # IMPORTANT: use kept names
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
    if MOS == "CCA":
        pycpt.plot_cca_modes(MOS, predictor_names, pxs, pys, domain_dir)
    pycpt.plot_forecasts(
        cpt_args, predictand_name, fcsts, domain_dir, predictor_names, MOS
    )

    # only attempt an mme if there is greater than 2 services picked up
    if len(predictor_names) < 2:
        warnings.warn(
            f"Skipping MME: only {len(predictor_names)} predictor available ({predictor_names}). "
            "MME requires at least 2 predictors."
        )
        print("Run complete. Outputs saved to:", case_dir)
        return

    ensemble = predictor_names

    det_fcst, pr_fcst, pev_fcst, nextgen_skill = pycpt.construct_mme(
        fcsts,
        hcsts,
        Y,
        ensemble,
        predictor_names,
        cpt_args,
        domain_dir,
    )

    mme_skill_metrics = [
        "spearman",
        "two_alternative_forced_choice",
        "root_mean_squared_error",
        "generalized_roc",
        "rank_probability_skill_score",
    ]

    pycpt.plot_mme_skill(
        predictor_names, nextgen_skill, MOS, domain_dir, mme_skill_metrics
    )
    pycpt.plot_mme_forecasts(
        cpt_args, predictand_name, pr_fcst, MOS, domain_dir, det_fcst
    )

    threshold = 0.5
    isPercentile = True

    (
        exceedance_prob,
        fcst_scale,
        climo_scale,
        fcst_mu,
        climo_mu,
        Y2,
        ntrain,
        transformed_threshold,
    ) = pycpt.construct_flex_fcst(
        MOS, cpt_args, det_fcst, threshold, isPercentile, Y, pev_fcst
    )

    if "station" in Y.coords:
        print(Y["Name"].to_pandas().to_string())

    # location_selector = {'X': 42.75, 'Y': -75.88}  # example for gridded data
    # location_selector = {'station': '11602'}  # example for station data
    location_selector = None
    forecast_colors = None

    pycpt.plot_mme_flex_forecast_v2(
        predictand_name,
        exceedance_prob,
        transformed_threshold,
        fcst_scale,
        climo_scale,
        fcst_mu,
        climo_mu,
        Y,
        Y2,
        ntrain,
        MOS,
        domain_dir,
        location_selector=location_selector,
        color_bar=forecast_colors,
    )

    print("Run complete. Outputs saved to:", case_dir)
    return
