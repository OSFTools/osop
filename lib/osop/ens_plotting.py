# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Collection of plotting codes relevant to ensembles."""

import logging

logger = logging.getLogger(__name__)

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from osop.plot_verify import location


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """Given a colormap, truncate it at bottom (minval) or top (maxval).

    Parameters
    ----------
    cmap : matplotlib.colors.Colormap
        Matplotlib color map to be truncated.
    minval : float, optional
        Minimum value to truncate at (0 to 1), by default 0.0.
    maxval : float, optional
        Maximum value to truncate at (0 to 1), by default 1.0.
    n : int, optional
        Number of color levels, by default 100.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        New truncated colormap.
    """
    new_cmap = colors.LinearSegmentedColormap.from_list(
        "trunc({n},{a:.2f},{b:.2f})".format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)),
    )

    return new_cmap


def get_cmap(precip_cs=False, wmo_cs=True):
    """Return 3 colormaps for below, normal and above tercile forecasts.

    Parameters
    ----------
    precip_cs : bool, optional
        If True, return a colorscale for rainfall, else for temperature, by default False.
    wmo_cs : bool, optional
        If True, use WMO colorscales, else use custom, by default True.

    Returns
    -------
    tuple
        Tuple of three colormaps: (below, normal, above).
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
    """Create a title string for the forecast plot based on the configuration.

    Parameters
    ----------
    config : dict
        Dictionary containing configuration parameters.

    Returns
    -------
    str
        Formatted title string.
    """
    lead = int(config["i"]) + 1
    atitle = (
        f"Tercile Summary:{config['origin']} {config['systemfc']} {config['fcstarty']} \n"
        f"startmonth:{config['start_month']} \nlead:{lead} \nVariable: {config['hc_var']}"
    )
    return atitle


def plot_tercile_fc(mme, atitle, var="precipitation", mask=None, map_setting="False"):
    """Plot a tercile forecast.

    Uses different colormaps for each of three terciles. Uses a threshold of 40% below which it does not plot.

    Parameters
    ----------
    mme : xarray.Dataset
        Tercile forecast dataset.
    atitle : str
        Title for the plot.
    var : str, optional
        Variable name in the dataset, by default "precipitation".
    mask : xarray.DataArray, optional
        Optional dry mask as a DataArray, by default None.
    map_setting : object or str, optional
        Optional map feature from cartopy, by default "False".

    Returns
    -------
    matplotlib.figure.Figure
        The matplotlib figure object.
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
    """Reformatt a forecast_percentage dataset to be able to run through ens_plotting routines.

    Parameters
    ----------
    data : xarray.Dataset
        The forecast percentage dataset.
    variable : str
        The variable name to use in the reformatted dataset.

    Returns
    -------
    xarray.Dataset
        A reformatted version of the forecast data for ens_plotting functions.
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
    except NameError as e:
        logger.error("Check data configurations match input parameters")
        raise e

    return new_dataset


def plot_forecasts(productdir, plotsdir, config):
    """Call functions and parse configurations for plotting forecasts.

    Parameters
    ----------
    productdir : str
        Location for dataset to be plotted.
    plotsdir : str
        Location for plots to save to.
    config : dict
        Dictionary for parameters of the file/dataset.

    Returns
    -------
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
        logger.info(f"Variable not identified: {variable}")
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
        plot_dataset_1m, atitle, var=variable, mask=None, map_setting=map_setting
    )
    # Save figure
    figname = f"{plotsdir}/{forecast_name_1m}.png"
    plt.savefig(figname, bbox_inches="tight", pad_inches=0.01)

    fig = plot_tercile_fc(
        plot_dataset_3m, atitle, var=variable, mask=None, map_setting=map_setting
    )
    # Save figure
    figname = f"{plotsdir}/{forecast_name_3m}.png"
    plt.savefig(figname, bbox_inches="tight", pad_inches=0.01)
