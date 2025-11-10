"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.

Collection of plotting codes relevant to ensembles
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib.colors as colors
import numpy as np
import xarray as xr
from osop.plot_verify import location


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """
    given a colormap, truncate it at bottom (minval)
    or top (maxval)

    Args:
        cmap - matplotlib color map to be truncated
        minval(float) - 0 to 1, min value to truncate at
        maxval(float) - 0 to 1, max value to truncate at

    Returns :
        new_cmap - new LinearSegmentedColormap
    """
    new_cmap = colors.LinearSegmentedColormap.from_list(
        "trunc({n},{a:.2f},{b:.2f})".format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)),
    )

    return new_cmap


def get_cmap(precip_cs=False, wmo_cs=True):
    """
    Returns 3 colormaps for below, normal and above
    tercile forecasts
    If precip_cs True, return a colorscale for rainfall,
    if false return colormap for temperature etc
    If wmo_cs true based on WMO colorscales used in
    LRF website else own version with single color
    """

    # use greys for normal tercile always
    # to have first color saturated use truncate_colormap
    cmap_normal = truncate_colormap(plt.cm.Greys, minval=0.2)
    if precip_cs:
        if wmo_cs:
            cmap_below = plt.cm.YlOrRd
            cmap_above = truncate_colormap(plt.cm.Greens, minval=0.2)
        else:
            # like WMO but with purples and browns
            cmap_below = truncate_colormap(
                plt.cm.BrBG, minval=0.0, maxval=0.15
            ).reversed()
            cmap_above = truncate_colormap(plt.cm.Purples, minval=0.5)
    else:
        # temperature etc
        if wmo_cs:
            # want start color a bit more saturatated for blues and greys
            cmap_below = truncate_colormap(plt.cm.Blues, minval=0.2)
            cmap_above = plt.cm.YlOrRd
        else:
            # like WMO but no yellow
            # want start color a bit more saturatated for blues and greys
            cmap_below = truncate_colormap(plt.cm.Blues, minval=0.2)
            cmap_above = truncate_colormap(plt.cm.Reds, minval=0.2)

    return cmap_below, cmap_normal, cmap_above

def fc_title(config):
    """
    Create a title string for the forecast plot based on the configuration.

    Parameters:
    - config: Dictionary containing configuration parameters.

    Returns:
    - atitle: Formatted title string.
    """
    lead = int(config["i"]) + 1
    atitle = (
        f"Tercile Summary:{config['origin']} {config['systemfc']} {config['fcstarty']} \n"
        f"startmonth:{config['start_month']} \nlead:{lead} \nVariable: {config['hc_var']}"
    )
    return atitle

def plot_tercile_fc(
    mme, atitle, var="precipitation", mask=None, map_setting="False"):
    """
    Function to plot a tercile forecast with different colormaps
    for each of three terciles. Uses a threshold of 40%
    below which it does not plot.

    Parameters:
    - mme: xarray DataArray with tercile forecasts
    - atitle: Title for the plot
    - plotsdir: Location to save the plot
    - var: Variable name in the dataset
    - mask: Optional dry mask as a DataArray
    - map_setting: Optional map feature from cartopy
    """

    # Apply threshold mask
    LTHRESH = 40.0
    Z1 = mme[var].data[0, ...]
    Z2 = mme[var].data[1, ...]
    Z3 = mme[var].data[2, ...]
    Z1_m = np.logical_or(np.logical_or(Z1 < Z2, Z1 < Z3), Z1 < LTHRESH)
    Z2_m = np.logical_or(np.logical_or(Z2 < Z1, Z2 < Z3), Z2 < LTHRESH)
    Z3_m = np.logical_or(np.logical_or(Z3 < Z1, Z3 < Z2), Z3 < LTHRESH)
    Z_all_m = np.stack([Z1_m, Z2_m, Z3_m])

    mme_mask = mme.copy(deep=True)
    mme_mask[var] = mme[var].where(np.logical_not(Z_all_m))

    # Choose colour scheme
    precip_cs = var in ["precip", "precipitation", "pr"]
    cmap_below, cmap_normal, cmap_above = get_cmap(precip_cs=precip_cs, wmo_cs=True)

    projection = ccrs.PlateCarree()
    fig = plt.figure()
    ax = plt.axes(projection=projection)

    levels = np.arange(LTHRESH, 85.0, 5.0)

    # Plot each tercile
    norm = BoundaryNorm(levels, ncolors=cmap_below.N, extend="max")
    clev1 = mme_mask[var][0, ...].plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap_below,
        norm=norm,
        add_colorbar=False,
    )

    norm = BoundaryNorm(levels, ncolors=plt.cm.Greys.N)
    clev2 = mme_mask[var][1, ...].plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap_normal,
        norm=norm,
        add_colorbar=False,
    )

    norm = BoundaryNorm(levels, ncolors=cmap_above.N, extend="max")
    clev3 = mme_mask[var][2, ...].plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap_above,
        norm=norm,
        add_colorbar=False,
    )

    # Set title once
    plt.title(atitle, fontsize=12, loc="left")
    ax.coastlines()

    # Optional map feature
    if map_setting != "False":
        ax.add_feature(map_setting, edgecolor="black", linewidth=0.5)

    # Colour bars
    cbar1 = plt.axes((0.15, 0.05, 0.2, 0.05))
    cbar = plt.colorbar(clev1, extend="max", location="bottom", cax=cbar1)
    cbar.set_label("Prob. below (%)")
    cbar.ax.invert_xaxis()
    for label in cbar.ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar2 = plt.axes((0.4, 0.05, 0.2, 0.05))
    cbar = plt.colorbar(clev2, location="bottom", cax=cbar2, extend="neither")
    cbar.set_label("Prob. near average (%)")
    for label in cbar.ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar3 = plt.axes((0.65, 0.05, 0.2, 0.05))
    cbar = plt.colorbar(clev3, extend="max", location="bottom", cax=cbar3)
    cbar.set_label("Prob. above (%)")
    for label in cbar.ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    # Gridlines
    ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=1,
        color="black",
        linestyle="--",
    )

    # Optional dry mask
    if mask is not None:
        mask.plot.contourf(
            ax=ax,
            transform=ccrs.PlateCarree(),
            add_colorbar=False,
            hatches=[".."],
            alpha=0,
            colors=None,
            add_labels=False,
        )
        plt.figtext(
            0.5,
            0.15,
            "Dotted area has average rainfall <0.1 mm/day",
            ha="center",
            fontsize=12,
        )

    return plt.gcf()


def reformatt(data, variable):
    """
    Takes a forecast_percentage dataset and reformats it to be able to run through
    ens_plotting routines.

    Parameters:
    config (dict): The cofiguraiton parameters for the forecast
    forecast_local (str): The location of the forecast_data set to be plotted

    Returns:
    new_dataset (xarray): A reformatted version of the forecast data for ens_plotting functions.
    """

    try:
        # Stack the three layers into a new dimension C
        values = np.stack(
            [data["lower"].values, data["middle"].values, data["higher"].values], axis=0
        )

        # Define new coordinates
        C = [1, 2, 3]  # 1: lower, 2: middle, 3: higher
        Y = data["lat"].values
        X = data["lon"].values

        # Create the new dataset
        new_dataset = xr.Dataset(
            {variable: (("C", "Y", "X"), values)}, coords={"C": C, "Y": Y, "X": X}
        )
    except NameError:
        print("Check data configurations match input parameters")

    return new_dataset


def plot_forecasts(productdir, plotsdir, config):
    """
    Calls functions and parses configurations over for plotting forecasts.

    Parameters:
    productdir (str): Location for dataset to be plot.
    plotsdir (str): Location for plots to save too.
    config (dict): Dictionary for parameters of the file/dataset.

    Returns:
    None

    """

    # forecast data set info
    forecast_local_1m = "{fpath}/{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.imonth_{i}.forecast_percentages.nc".format(
        fpath=productdir, **config
    )
    forecast_name_1m = "{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.imonth_{i}".format(
        **config
    )
    forecast_local_3m = "{fpath}/{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.3m.forecast_percentages.nc".format(
        fpath=productdir, **config
    )
    forecast_name_3m = "{origin}_{systemfc}_{fcstarty}-{fcendy}_monthly_mean_{start_month}_{leads_str}_{area_str}_{hc_var}.3m".format(
        **config
    )

    # checks the variable for use.
    variable = "{hc_var}".format(**config)
    if variable == "2m_temperature":
        variable = "temperature"
    elif variable == "total_precipitation":
        variable = "precipitation"
    else:
        print("Variable not identified")
        # Future functionality should be able to handle this - see plot_tercile_fc
        variable = variable

    # Reformat dataset for plotting
    fcst_local_1m = xr.open_dataset(forecast_local_1m)
    fcst_local_3m = xr.open_dataset(forecast_local_3m)
    fcst_local_3m = fcst_local_3m.squeeze(dim="forecastMonth")

    plot_dataset_1m = reformatt(fcst_local_1m, variable)
    plot_dataset_3m = reformatt(fcst_local_3m, variable)

    # Tercile Summary - 1month forecasts, per origin centre.
    map_setting = location(config)
    atitle = fc_title(config)
    fig = plot_tercile_fc(
        plot_dataset_1m, atitle, var=variable, mask=None,
        map_setting=map_setting
    )
        # Save figure
    figname = f"{plotsdir}/{forecast_name_1m}.png"
    plt.savefig(figname, bbox_inches="tight", pad_inches=0.01)

    fig = plot_tercile_fc(
        plot_dataset_3m,  atitle, var=variable, mask=None,
        map_setting=map_setting
    )
    # Save figure
    figname = f"{plotsdir}/{forecast_name_3m}.png"
    plt.savefig(figname, bbox_inches="tight", pad_inches=0.01)
