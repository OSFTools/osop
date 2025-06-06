# Libraries for working with multi-dimensional arrays
import xarray as xr
import numpy as np
import os

# Libraries for plotting and geospatial data visualisation
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy import crs as ccrs
import calendar



CATNAMES = ["lower tercile", "middle tercile", "upper tercile"]

BORDER_OPT={'Morocco':'admin_0_countries_mar',
            'UK':'admin_0_countries_gbr',
           'SAU':'admin_0_countries_sau'
           }

def location(config):
    """
    Prepares location specific POV borders for plots based on arguments in the config dictonary
    Args:
        config (dict): A dictionary containing the configuration parameters.
    Returns:
        A Natural Earth data set name to go into the axis plot that's downloaded based on location; if no location set 
        i.e. None - no borders will plot.
        If the name is misspelt then a key error will raise suggesting a check of locaiton entry in the shell script. 
    Redundancy:
        Natural Earth has a download issue that searches for a file that dosnt exist; a partial import is managed regardless
        On second run - as this file is not used and download has already happened -  the plot will work fine.
        To avoid a second run each time a new data set is imported the try/except does the import for no reason and then the finally is used after 
        to generate the plot. - This is a Natural Earth Specific problem that can be removed when fixed.
        Relevent to Cartopy issue #2319 , #2477 amd #2534 - when resolved can be removed 
    """
    if config['border'] in BORDER_OPT:
        border_set = BORDER_OPT[config['border']]
        try:
            shpfilename = shpreader.natural_earth(resolution='10m', category='cultural', name=border_set)
        except KeyError:
            shpfilename = shpreader.natural_earth(resolution='10m', category='cultural', name=border_set)
        finally:
            local = cfeature.NaturalEarthFeature(
                category='cultural',
                name=border_set,
                scale='10m',
                facecolor='none')
            return local
    elif config['border'] == 'None':
        local = 'False'
        return local
    else:
        raise KeyError("Location Name does not exist in dictionary. Please check spelling of location input or type None for no borders.")
    
        
    


def prep_titles(config):
    """
    Prepare titles for the plot based on the arguments in the config dictionary.
    Currently this replicates the titles here:
    https://confluence.ecmwf.int/display/CKB/C3S+seasonal+forecasts+verification+plots

    Args:
        config (dict): A dictionary containing the configuration parameters.

    Returns:
        tuple: A tuple containing the prepared titles for the plot.
            The first element is the first line of the title.
            The second element is the second line of the title.
            The third elementis the third line of the title.
    """
    tit_line1 = "{origin} {system}".format(**config)
    tit_line2_base = (
        f'Start month: {calendar.month_abbr[config["start_month"]].upper()}'
    )

    if config["aggr"] == "1m":
        validmonth = config["valid_month"]
        validmonth = validmonth if validmonth <= 12 else validmonth - 12
        tit_line2 = (
            tit_line2_base
            + f" - Valid month: {calendar.month_abbr[validmonth].upper()}"
        )
    elif config["aggr"] == "3m":
        validmonths = [
            vm if vm <= 12 else vm - 12
            for vm in [config["valid_month"] + shift for shift in range(3)]
        ]
        validmonths = [calendar.month_abbr[vm][0] for vm in validmonths]
        tit_line2 = tit_line2_base + f' - Valid months: {"".join(validmonths)}'
    else:
        raise BaseException(f"Unexpected aggregation")
    tit_line3 = f"Verification period: {config['hcstarty']} - {config['hcstarty']}"
    return tit_line1, tit_line2, tit_line3


def plot_score(score_f, score_fname, category, config, score, titles, datadir):
    """
    Plot the score on a map.

    Parameters:
    score_f (numpy.ndarray): The score data.
    category (str): The tercile category. Set to None if not relevant for the score.
    config (dict): Configuration parameters.
    score (str): The name of the score.
    titles (list): Titles for the plot.
    datadir (str): The directory to save the plot.

    Returns:
    None
    """
    # choose a projection and set up an axes based on that
    projection = ccrs.Mercator()
    ax = plt.axes(projection=projection)
    # Specify CRS
    crs = ccrs.PlateCarree()

    info = f'{score_fname}_category_{category}'
    if score == "rps":
        p = score_f[config["var"]][0, :, :]
        lon = score_f[config["var"]].lon
        lat = score_f[config["var"]].lat
        cols = "YlGn_r"
        ex_dir = "max"
        under = "purple"
        levels = np.linspace(0.0, 0.5, 11)
        info = score_fname
        plt.title(
            f"{score} \n" + titles[0] + f" {config['var']}\n" + 
            titles[1] + "\n" + titles[2], loc="left"
        )
    elif score == "bs":
        p = score_f[config["var"]].sel(category=category)[0, :, :]
        lon = score_f[config["var"]].sel(category=category).lon
        lat = score_f[config["var"]].sel(category=category).lat
        cols = "YlGn_r"
        ex_dir = "max"
        under = "purple"
        plt.title(
            f"{score} \n"
            + titles[0]
            + f' {config["var"]}'
            + f" ({CATNAMES[category]})\n"
            + titles[1] + "\n" + titles[2],
            loc="left",
        )
        levels = np.linspace(0.0, 0.5, 11)
    elif score in ["roc", "rocss"]:
        plt.title(
            f"{score} \n"
            + titles[0]
            + f' {config["var"]}'
            + f" ({CATNAMES[category]})\n"
            + titles[1] + "\n" + titles[2],
            loc="left",
        )
        p = score_f[config["var"]].sel(category=category)[0, :, :]
        lon = score_f[config["var"]].sel(category=category).lon
        lat = score_f[config["var"]].sel(category=category).lat
        cols = "YlGn"
        ex_dir = "neither"
        under = "lightgray"
        levels = np.linspace(0.5, 1.0, 6)
    
    info = f'{score_fname}_category_{category}'
    if score == 'rps':
        p = score_f[config['var']][0,:,:]
        lon = score_f[config['var']].lon
        lat = score_f[config['var']].lat
        cols = 'YlGn_r'
        ex_dir = 'max'
        under = 'purple'
        levels = np.linspace(0.,0.5,11)
        info = score_fname
        plt.title(f'{score} \n' + titles[0] + f' {config["var"]}\n' 
                  + titles[1] + '\n' + titles[2], loc='left')
    elif score == 'bs':
        p = score_f[config['var']].sel(category=category)[0,:,:]
        lon = score_f[config['var']].sel(category=category).lon
        lat = score_f[config['var']].sel(category=category).lat
        cols = 'YlGn_r'
        ex_dir = 'max'
        under = 'purple'
        plt.title(f'{score} \n' + titles[0] + f' {config["var"]}' 
                  + f' ({CATNAMES[category]})\n' + titles[1]
                  + '\n' + titles[2],  loc='left')
        levels = np.linspace(0.,0.5,11)
    elif score in ['roc', 'rocss']:
        plt.title(f'{score} \n' + titles[0] + f' {config["var"]}' + 
                  f' ({CATNAMES[category]})\n' + titles[1] +
                  '\n' + titles[2],  loc='left')
        p = score_f[config['var']].sel(category=category)[0,:,:]
        lon = score_f[config['var']].sel(category=category).lon
        lat = score_f[config['var']].sel(category=category).lat
        cols = 'YlGn'
        ex_dir = 'neither'
        under = 'lightgray'
        levels = np.linspace(0.5,1.,6)

    cs = plt.contourf(
        lon,
        lat,
        p,
        transform=ccrs.PlateCarree(),
        levels=levels,
        cmap=cols,
        extend=ex_dir,
    )
    
    
    cs.cmap.set_under(under)

    map_setting = location(config)
    if map_setting != "False":
        ax.add_feature(map_setting, edgecolor="black", linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, edgecolor="black", linewidth=2.0)


    print(info)
    plt.colorbar()
    plt.savefig(os.path.join(datadir, "scores", f"{info}.png"))
    plt.close()


def plot_rel(score_f, score_fname, config, score, datadir, titles):
    """Plot reliability diagram
    Parameters:
    score_f (numpy.ndarray): The reliability score data.
    config (dict): Configuration parameters.
    score (str): The name of the score.
    datadir (str): The directory to save the plot.

    Returns:
    None
    """
    info = score_fname

    for var in score_f.data_vars:
        p = score_f[var][0, :, :]

        fig, ax1 = plt.subplots(figsize=(18, 10))
        # ax2 = ax1.twinx()

        ax1.plot(
            [0, 1], [0, 1], linestyle="--", color="gray", label="Perfect Reliability"
        )
        for icat in p.category.values:
            cat_lab = f" {CATNAMES[icat]}"
            avalues = p.sel(category=icat)
            ax1.plot(avalues.forecast_probability.values, avalues.values, label=cat_lab)
            # ax2.hist(avalues.values, bins=10, alpha=0.5, label=f'category {icat}')

        ax1.set_xlabel("Forecast Probability")
        ax1.set_ylabel("Values")
        # ax2.set_ylabel('Histogram')

        ax1.legend(loc="upper left")
        # ax2.legend(loc='upper right')
        plt.title(
            titles[0] + f" {config['var']}\n" + 
            titles[1] + "\n" + titles[2], loc="left"
        )
        plt.tight_layout()
        plt.savefig(os.path.join(datadir, "scores", f"{score_fname}.png"))
        plt.close()


def corr_plots(datadir, hcst_bname, aggr, config, score, titles):
    """Plot deterministic scores
    Parameters:
    datadir (str): The directory to save the plot.
    hcst_bname (str): The basename of the hindcast file.
    aggr (str): The aggregation period.
    config (dict): Configuration parameters.
    titles (list): Titles for the plot.

    Returns:
    None
    """
    # Read the data files
    corr = xr.open_dataset(f"{datadir}/scores/{hcst_bname}.{aggr}.{score}.nc")
    corr_pval = xr.open_dataset(f"{datadir}/scores/{hcst_bname}.{aggr}.{score}_pval.nc")

    # Rearrange the dataset longitude values for plotting purposes
    corr = corr.assign_coords(lon=(((corr.lon + 180) % 360) - 180)).sortby("lon")
    corr_pval = corr_pval.assign_coords(
        lon=(((corr_pval.lon + 180) % 360) - 180)
    ).sortby("lon")

    # Set up the figure and axes
    
    fig = plt.figure(figsize=(18, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    map_setting = location(config)
    if map_setting != "False":
        ax.add_feature(map_setting, edgecolor="black", linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, edgecolor="black", linewidth=2.0)

        

        

    corrvalues = corr[config["var"]][0, :, :].values
    corrpvalvalues = corr_pval[config["var"]][0, :, :].values

    if corrvalues.T.shape == (
        corr[config["var"]].lat.size,
        corr[config["var"]].lon.size,
    ):
        print("Data values matrices need to be transposed")
        corrvalues = corrvalues.T
        corrpvalvalues = corrpvalvalues.T
    elif corrvalues.shape == (
        corr[config["var"]].lat.size,
        corr[config["var"]].lon.size,
    ):
        pass
    else:
        raise BaseException(f"Unexpected data value matrix shape: {corrvalues.shape}")
    plt.contourf(
        corr[config["var"]].lon,
        corr[config["var"]].lat,
        corrvalues,
        levels=np.linspace(-1.0, 1.0, 11),
        cmap="RdYlBu_r",
    )
    cb = plt.colorbar(shrink=0.5)
    cb.ax.set_ylabel({score}, fontsize=12)
    origylim = ax.get_ylim()

    # add stippling for significance 
    # where p value > 0.05
    # note can get NaN values in the pval matrix
    # where the standard deviation of one of the 
    # fields is zero. Use nanmax not max
    
    plt.contourf(
        corr_pval[config["var"]].lon,
        corr_pval[config["var"]].lat,
        corrpvalvalues,
        levels=[np.nanmin(corrpvalvalues), 0.05, np.nanmax(corrpvalvalues)],
        hatches=["", "..."],
        colors="none",
    )
    # We need to ensure after running plt.contourf() the ylim hasn't changed
    if ax.get_ylim() != origylim:
        ax.set_ylim(origylim)

    plt.title(
        f"{titles[0]} (stippling where p>0.05)\n {titles[1]} \n {titles[2]}",
        loc="left",
    )
    plt.tight_layout()
    figname = f"{datadir}/scores/{hcst_bname}.{aggr}.{score}.png"
    plt.savefig(figname)
    plt.close()



def generate_plots(config, titles, downloaddir):
    ## read in the data
    score_fname = "{origin}_{system}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{fname_var}.{aggr}.{score}.nc".format(
        **config
    )
    score_data = xr.open_dataset(os.path.join(downloaddir, "scores", score_fname))

    if config["score"] == "spearman_corr":
        score_fname = "{origin}_{system}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{fname_var}".format(
            **config
        )
        corr_plots(downloaddir, score_fname, config["aggr"], config, config["score"], titles)

    elif config["score"] == "pearson_corr":
        score_fname = "{origin}_{system}_{hcstarty}-{hcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{fname_var}".format(
            **config
        )
        corr_plots(downloaddir, score_fname, config["aggr"], config, config["score"], titles)

    elif config["score"] == "rel":
        plot_rel(score_data, score_fname, config, config["score"], downloaddir, titles)

    elif config["score"] == "rps":
        plot_score(score_data, score_fname, None, config, config["score"], titles, downloaddir)

    elif config["score"] == "bs":
        for cat in [0, 1, 2]:
            plot_score(score_data, score_fname, cat, config, config["score"], titles, downloaddir)

    else:
        for cat in [0, 1, 2]:
            config["category"] = cat
            plot_score(score_data, score_fname, cat, config, config["score"], titles, downloaddir)
