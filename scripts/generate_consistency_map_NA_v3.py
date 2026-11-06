#!/usr/bin/env python
"""
WMO LC-LRFMME Consistency Map — JJA 2026 Precipitation — v2
============================================================
Adaptation aux noms de coordonnées GRIB (g0_lat_1, g0_lon_2)
"""
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

# ============================================================================
NC_PATH = r"C:/Users/wafae/Desktop/DATA PC HP/RCOFS/MEDCOF/MEDCOF26 JJA2026/cons_fcst_202605_202606_202608_JJA_prec.nc"
OUT_DIR = "C:/Users/wafae/projets/osop-data/seafoam/data/master/RCC_NA_maps"
os.makedirs(OUT_DIR, exist_ok=True)

EXTENT = [-18, 38, 14, 38]  # Afrique du Nord

WMO_CENTRES = "Beijing, CMCC, CPTEC, ECMWF, Exeter, Melbourne, Montreal, " \
              "Moscow, Offenbach, Pune, Seoul, Tokyo, Toulouse, Washington"

# ============================================================================
# FRONTIÈRES MAROC
# ============================================================================
def get_morocco_borders():
    border_set = "admin_0_countries_mar"
    try:
        _ = shpreader.natural_earth(resolution="10m", category="cultural", name=border_set)
    except KeyError:
        _ = shpreader.natural_earth(resolution="10m", category="cultural", name=border_set)
    finally:
        local = cfeature.NaturalEarthFeature(
            category="cultural", name=border_set, scale="10m", facecolor="none")
    return local

print("Chargement shapefile Maroc...")
MA_BORDERS = get_morocco_borders()
print("✅\n")

# ============================================================================
# PALETTE WMO
# ============================================================================
WMO_COLORS = [
    "#08306B",  # -14 dark blue
    "#1F5BA5",  # -13
    "#4B8FD0",  # -12
    "#9CCAE8",  # -11
    "#D6EBF5",  # -10 light blue
    "#FFFFFF",  # gap zone (-9 to +9) : NO CONSENSUS = white
    "#FFFACC",  # +10 pale yellow
    "#FED98E",  # +11
    "#FB9034",  # +12
    "#E0451E",  # +13
    "#A50026",  # +14 dark red
]
WMO_CMAP = ListedColormap(WMO_COLORS)
WMO_BOUNDS = [-14.5, -13.5, -12.5, -11.5, -10.5, -9.5,
              9.5, 10.5, 11.5, 12.5, 13.5, 14.5]
WMO_NORM = BoundaryNorm(WMO_BOUNDS, ncolors=len(WMO_COLORS), extend='neither')

# ============================================================================
# 1. CHARGEMENT
# ============================================================================
print("=" * 70)
print("LECTURE DU NETCDF")
print("=" * 70)

ds = xr.open_dataset(NC_PATH)
print(f"Variables : {list(ds.data_vars)}")
print(f"Coords   : {list(ds.coords)}")

# Noms GRIB-style spécifiques
LAT_NAME = 'g0_lat_1'
LON_NAME = 'g0_lon_2'
VAR_NAME = 'prec'

# Champ
field = ds[VAR_NAME]
print(f"\nField stats : min={float(field.min()):.1f}, max={float(field.max()):.1f}")

# Convertir longitudes 0-360 → -180/180
lon_vals = ds[LON_NAME].values
if lon_vals.max() > 180:
    print("\nConversion longitudes 0–360 → -180/180...")
    new_lon = np.where(lon_vals > 180, lon_vals - 360, lon_vals)
    ds = ds.assign_coords({LON_NAME: new_lon}).sortby(LON_NAME)
    field = ds[VAR_NAME]
    print("✅ Conversion faite")

lat = ds[LAT_NAME].values
lon = ds[LON_NAME].values
print(f"\nFinal lat range : {lat.min():.1f}..{lat.max():.1f}°N")
print(f"Final lon range : {lon.min():.1f}..{lon.max():.1f}°E")

# ============================================================================
# 2. CARTE
# ============================================================================
print("\n" + "=" * 70)
print("GÉNÉRATION CARTE")
print("=" * 70)

fig = plt.figure(figsize=(13, 7))

# Carte
ax = fig.add_axes([0.05, 0.18, 0.90, 0.68], projection=ccrs.PlateCarree())
ax.set_extent(EXTENT, crs=ccrs.PlateCarree())

# Données
im = ax.pcolormesh(lon, lat, field.values,
                    transform=ccrs.PlateCarree(),
                    cmap=WMO_CMAP, norm=WMO_NORM, shading='auto',
                    zorder=2)

# Océan PAR-DESSUS
ax.add_feature(cfeature.OCEAN, facecolor='#E8F4F8', zorder=10)
ax.add_feature(cfeature.LAKES, facecolor='#E8F4F8', zorder=10)
ax.add_feature(cfeature.COASTLINE, linewidth=0.7, edgecolor='black', zorder=11)
ax.add_feature(MA_BORDERS, linewidth=0.5, edgecolor='black', zorder=11)

# Grille
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray',
                   linestyle='--', alpha=0.5, zorder=12)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 9}
gl.ylabel_style = {'size': 9}

# === TITRES STYLE WMO ===
fig.text(0.05, 0.93, 'Consistency Map',
          ha='left', fontsize=14, fontweight='bold')
fig.text(0.05, 0.905, WMO_CENTRES,
          ha='left', fontsize=8, style='italic', color='#444444')
fig.text(0.05, 0.875, 'Precipitation : JJA 2026',
          ha='left', fontsize=12, fontweight='bold')
fig.text(0.95, 0.875, '(issued on May 2026)',
          ha='right', fontsize=10, style='italic')

# Colorbar
cax = fig.add_axes([0.15, 0.08, 0.70, 0.025])
cb = plt.colorbar(im, cax=cax, orientation='horizontal',
                   ticks=[-14, -13, -12, -11, -10, 10, 11, 12, 13, 14])
cb.ax.tick_params(labelsize=9)
# Label sous la colorbar pour expliquer le gap blanc
cb.set_label("Number of models agreeing on anomaly sign (white = no consensus)",
              fontsize=8, style='italic', labelpad=4)

# Footer
fig.text(0.5, 0.025,
         '** where the positive numbers mean the number of models that predict positive anomaly and vice versa **',
         ha='center', fontsize=8, style='italic', color='#333333')

fig.text(0.93, 0.06, 'WMO Lead Centre for\nLR-MME',
          ha='center', fontsize=7, color='#2C5AA0', style='italic')

out_file = os.path.join(OUT_DIR, "Consistency_Map_TP_JJA2026_NorthAfrica.png")
plt.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.2)
plt.close()

print(f"\n✅ {out_file}")
print(f"   {os.path.getsize(out_file)/1024:.0f} KB")

# Statistiques de signal sur la zone Afrique du Nord
print("\n" + "=" * 70)
print("STATISTIQUES SUR L'AFRIQUE DU NORD")
print("=" * 70)
field_na = field.sel({LAT_NAME: slice(14, 38), LON_NAME: slice(-18, 38)})
n_total = field_na.size
n_strong_neg = int((field_na <= -10).sum())   # consensus sec
n_strong_pos = int((field_na >= 10).sum())    # consensus humide
n_no_signal = n_total - n_strong_neg - n_strong_pos
print(f"  Total pixels NA       : {n_total}")
print(f"  Consensus sec  (≤-10) : {n_strong_neg} ({100*n_strong_neg/n_total:.1f}%)")
print(f"  Consensus humide (≥10): {n_strong_pos} ({100*n_strong_pos/n_total:.1f}%)")
print(f"  Pas de consensus      : {n_no_signal} ({100*n_no_signal/n_total:.1f}%)")

ds.close()
print("\n" + "=" * 70)
print(f"  explorer \"{OUT_DIR}\"")
