#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Carte de consistance LC-LRFMME recadree sur le domaine RCC Afrique du Nord.
Valeurs = +/-(N + 0.1), N = nombre de modeles d'accord sur le signe de l'anomalie.
"""
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

SRC = r"C:\Users\wafae\Downloads\cons_fcst_202608_202609_202611_SON_prec.nc"
DST = (r"C:\Users\wafae\projets\osop-data\seafoam\data\master\RCC_NA_maps"
       r"\Consistency_Map_TP_SON2026_NorthAfrica.png")
EXTENT = [-20, 40, 12, 40]

ds = xr.open_dataset(SRC)
da = ds["prec"].rename({"g0_lat_1": "lat", "g0_lon_2": "lon"})
da = da.assign_coords(lon=(((da.lon + 180) % 360) - 180)).sortby("lon")
da = da.sel(lat=slice(EXTENT[2] - 3, EXTENT[3] + 3),
            lon=slice(EXTENT[0] - 3, EXTENT[1] + 3))

cols = ["#08306b", "#2171b5", "#4292c6", "#6baed6", "#c6dbef",
        "#ffffff",
        "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02"]
bounds = [-14.5, -13.5, -12.5, -11.5, -10.5, -9.5,
          9.5, 10.5, 11.5, 12.5, 13.5, 14.5]
cmap, norm = ListedColormap(cols), BoundaryNorm(bounds, 11)

fig = plt.figure(figsize=(11, 5.2))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
m = ax.pcolormesh(da.lon, da.lat, da.values, cmap=cmap, norm=norm,
                  transform=ccrs.PlateCarree(), shading="nearest")
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linewidth=0.6, edgecolor="0.25")
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="gray", alpha=0.4)
gl.top_labels = gl.right_labels = False
gl.xlabel_style = gl.ylabel_style = {"size": 9}

cb = fig.colorbar(m, ax=ax, orientation="horizontal", pad=0.09, shrink=0.8,
                  ticks=[-14, -13, -12, -11, -10, 10, 11, 12, 13, 14])
cb.set_label("Number of models agreeing on the sign of the anomaly", fontsize=9)
ax.set_title("Consistency Map - Precipitation SON 2026\n"
             "WMO LC-LRFMME, issued August 2026 - 14 contributing centres",
             fontsize=12, fontweight="bold", pad=12)

fig.savefig(DST, dpi=200, bbox_inches="tight", facecolor="white")
print("[OK] " + DST)

sub = da.sel(lat=slice(14, 38), lon=slice(-18, 37)).values
print(f"  humide >= 10 centres : {(sub >= 9.5).mean()*100:.0f} % du domaine")
print(f"  sec    <= -10        : {(sub <= -9.5).mean()*100:.0f} %")
