"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.
"""

# Ensure the top level directory has been added to PYTHONPATH

# Import functions 
import os
import yaml
from yaml.loader import SafeLoader

#import needed local functions
from osop.plot_verify import generate_plots, prep_titles


if __name__ == "__main__":
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the plotting functions to generate verification plots
    """
    scores = ["spearman_corr","pearson_corr", "roc", "rocss", "rps", "rel", "bs"]

    # get command line args
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
            location_test = config_test["location"]
            downloaddir_test = config_test["downloaddir"]
            scoresdir = config_test["scoresdir"]
            plotsdir = config_test["plotdir"]
            years = config_test["years"]
            centre = config_test["centre"]
            print('yml success')
        except yaml.YAMLError as e:
            print(e)
    
    downloaddir_test = os.path.expandvars(downloaddir_test)
    scoresdir = os.path.expandvars(scoresdir)
    plotsdir = os.path.expandvars(plotsdir)
    os.makedirs(plotsdir, exist_ok=True)

    # unpack args and reformat if needed
    border = location_test
    downloaddir = downloaddir_test
    scoresdir = scoresdir
    plotdir = plotsdir
    month = int(Month_test)
    leads = leads_test
    leadtime_month = [int(l) for l in leads_test.split(",")]
    leads_str = "".join([str(mon) for mon in leadtime_month])
    obs_str = "".join([str(mon - 1) for mon in leadtime_month])
    area = [float(pt) for pt in area_test.split(",")]
    area_str = area_test.replace(",", ":")
    fname_var = varaible_test

    if fname_var == "2m_temperature":
        var = "t2m"
    elif fname_var == "total_precipitation":
        var = "tprate"
    else:
        raise ValueError(f"Unknown hindcast variable: {fname_var}")

    valid_month = month + (leadtime_month[0] - 1)

    # add arguments to config
    config = dict(
        border = border,
        start_month=month,
        valid_month=valid_month,
        area_str=area_str,
        leads_str=leads_str,
        obs_str=obs_str,
        fname_var=fname_var,
        var=var,
    )


    if years:
        config["hcstarty"] = years[0]
        config["hcendy"] = years[1]
    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016

    for centre in config_test["centre"]:
        config["origin"] = centre
        for score in scores:
            for aggr in ["1m", "3m"]:
                config["aggr"] = aggr
                config["score"] = score

                # run for appropriate system
                if centre == "eccc":
                    # two models aka systems are live - call twice with each system number
                    config["system"] = Services["eccc_can"]
                    ## set titles
                    titles = prep_titles(config)
                    generate_plots(config, titles, scoresdir, plotdir)

                    ## repeat for second system
                    config["system"] = Services["eccc_gem5"]
                    ## set titles
                    titles = prep_titles(config)
                    generate_plots(config, titles, scoresdir, plotdir)
                else:
                    if centre not in Services.keys():
                        raise ValueError(f"Unknown system for C3S: {centre}")
                    config["system"] = Services[centre]
                    ## set titles
                    titles = prep_titles(config)
                    generate_plots(config, titles, scoresdir, plotdir)
