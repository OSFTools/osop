# Ensure the top level directory has been added to PYTHONPATH
import argparse
#import needed local functions
from osop.constants import SYSTEMS
from osop.plot_verify import generate_plots, prep_titles




def parse_args():
    """
    set up argparse to get command line arguments

    Returns:
        args: argparse args object
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--location", required=True)
    parser.add_argument("--centre", required=True, help="centre to download")
    parser.add_argument("--month", required=True, help="start month for hindcasts")
    parser.add_argument(
        "--variable",
        required=True,
        help="variable to verify. 2m_temperature or total_precipitation",
    )
    parser.add_argument(
        "--leads", required=True, help="forecast range in months (comma separated)"
    )
    parser.add_argument(
        "--area",
        required=True,
        help="sub-area in degrees for retrieval (comma separated N,W,S,E)",
    )
    parser.add_argument("--downloaddir", required=True, help="location to download to")
    parser.add_argument(
        "--years",
        required=False,
        help="Years to rerieve data for (comma separated). Optional. Default is hindcast period 1993-2016.",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    """
    Called when this is run as a script
    Get the command line arguments using argparse
    Call the plotting functions to generate verification plots
    """
    scores = ["spearman_corr", "roc", "rocss", "rps", "rel", "bs"]

    # get command line args
    args = parse_args()

    # unpack args and reformat if needed
    border = args.location
    centre = args.centre
    downloaddir = args.downloaddir
    month = int(args.month)
    leads = args.leads
    leadtime_month = [int(l) for l in args.leads.split(",")]
    leads_str = "".join([str(mon) for mon in leadtime_month])
    obs_str = "".join([str(mon - 1) for mon in leadtime_month])
    area = [float(pt) for pt in args.area.split(",")]
    area_str = args.area.replace(",", ":")
    fname_var = args.variable

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
        origin=centre,
        area_str=area_str,
        leads_str=leads_str,
        obs_str=obs_str,
        fname_var=fname_var,
        var=var,
    )

    if args.years:
        config["hcstarty"] = args.years[0]
        config["hcendy"] = args.years[1]
    else:
        config["hcstarty"] = 1993
        config["hcendy"] = 2016

    for score in scores:
        for aggr in ["1m", "3m"]:
            config["aggr"] = aggr
            config["score"] = score

            # run for appropriate system
            if centre == "eccc":
                # two models aka systems are live - call twice with each system number
                config["system"] = SYSTEMS["eccc_can"]
                ## set titles
                titles = prep_titles(config)
                generate_plots(config, titles, downloaddir)

                ## repeat for second system
                config["system"] = SYSTEMS["eccc_gem5"]
                ## set titles
                titles = prep_titles(config)
                generate_plots(config, titles, downloaddir)
            else:
                if centre not in SYSTEMS.keys():
                    raise ValueError(f"Unknown system for C3S: {centre}")
                config["system"] = SYSTEMS[centre]
                ## set titles
                titles = prep_titles(config)
                generate_plots(config, titles, downloaddir)
