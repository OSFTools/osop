#!/usr/bin/env python
"""
OSOP/SEAFOAM Multi-System TP Forecast — SON 2026 — North Africa
================================================================
v3 : Numéros de système mis à jour (aout 2026)
  - ukmo 604 → 610
  - cmcc 35 → 4 (SPS 3.5 → SPS 4.0)
  - jma 3 → 4
"""
import os
import glob
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm
from matplotlib.gridspec import GridSpec
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

FCST_DIR = "C:/Users/wafae/projets/osop-data/seafoam/data/master/forecast/products"
OUT_DIR = "C:/Users/wafae/projets/osop-data/seafoam/data/master/RCC_NA_maps"
os.makedirs(OUT_DIR, exist_ok=True)

EXTENT = [-18, 38, 14, 38]
THRESHOLD = 40.0

# ============================================================================
# MAPPING MIS À JOUR (aout 2026)
# ============================================================================
MODEL_LABELS = {
    'ecmwf_51':       'ECMWF SEAS5',
    'ukmo_610':       'UK Met Office GloSea6.1',
    'meteo_france_9': 'Meteo-France System 9',
    'dwd_22':         'DWD GCFS 2.1',
    'cmcc_4':         'CMCC SPS 4.0',
    'ncep_2':         'NCEP CFSv2',
    'eccc_4':         'ECCC CanSIPSv3',
    'eccc_5':         'ECCC GEM5',
    'jma_4':          'JMA CPS4 (excl. MME)',
    'mme_1':          'MULTI-MODEL ENSEMBLE',
}
MODEL_ORDER = ['ecmwf_51', 'ukmo_610', 'meteo_france_9', 'dwd_22',  'cmcc_4',
               'ncep_2',   'eccc_4',   'eccc_5',          'jma_4',   'mme_1']

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
print("OK\n")

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    return mcolors.LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n))
    )

TP_CMAP_BELOW  = truncate_colormap(plt.cm.YlOrBr, minval=0.3)
TP_CMAP_NORMAL = truncate_colormap(plt.cm.Greys,  minval=0.2, maxval=0.5)
TP_CMAP_ABOVE  = truncate_colormap(plt.cm.Greens, minval=0.3)

LEVELS = np.arange(THRESHOLD, 85.0, 5.0)
NORM_BELOW  = BoundaryNorm(LEVELS, ncolors=TP_CMAP_BELOW.N,  extend='max')
NORM_NORMAL = BoundaryNorm(LEVELS, ncolors=TP_CMAP_NORMAL.N, extend='max')
NORM_ABOVE  = BoundaryNorm(LEVELS, ncolors=TP_CMAP_ABOVE.N,  extend='max')

def find_file(prefix):
    pattern = f"{prefix}_2026-2026*monthly_mean_8_234_45_-30_-2.5_60_total_precipitation*3m*forecast_percentages.nc"
    matches = glob.glob(os.path.join(FCST_DIR, pattern))
    return matches[0] if matches else None

def setup_subplot(ax, title=''):
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, facecolor='#E8F4F8', zorder=0)
    ax.add_feature(cfeature.LAND,  facecolor='#FAFAF5', zorder=0)
    if title:
        ax.set_title(title, fontsize=10, fontweight='bold', pad=4)

def plot_tercile(ax, ds):
    ds = ds.isel(forecastMonth=0)
    lower  = ds['lower'].values
    middle = ds['middle'].values
    higher = ds['higher'].values

    mask_b = np.logical_or(np.logical_or(lower  < middle, lower  < higher), lower  < THRESHOLD)
    mask_n = np.logical_or(np.logical_or(middle < lower,  middle < higher), middle < THRESHOLD)
    mask_a = np.logical_or(np.logical_or(higher < lower,  higher < middle), higher < THRESHOLD)

    fb = np.where(~mask_b, lower,  np.nan)
    fn = np.where(~mask_n, middle, np.nan)
    fa = np.where(~mask_a, higher, np.nan)

    im_b = ax.pcolormesh(ds.lon, ds.lat, fb, transform=ccrs.PlateCarree(),
                          cmap=TP_CMAP_BELOW, norm=NORM_BELOW, shading='auto', zorder=2)
    im_n = ax.pcolormesh(ds.lon, ds.lat, fn, transform=ccrs.PlateCarree(),
                          cmap=TP_CMAP_NORMAL, norm=NORM_NORMAL, shading='auto', zorder=2)
    im_a = ax.pcolormesh(ds.lon, ds.lat, fa, transform=ccrs.PlateCarree(),
                          cmap=TP_CMAP_ABOVE, norm=NORM_ABOVE, shading='auto', zorder=2)

    ax.add_feature(cfeature.OCEAN, facecolor='#E8F4F8', zorder=10)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor='black', zorder=11)
    ax.add_feature(MA_BORDERS, linewidth=0.4, edgecolor='black', zorder=11)
    return im_b, im_n, im_a

print("=" * 70)
print("RECHERCHE DES FICHIERS OSOP (versions aout 2026)")
print("=" * 70)
datasets = {}
for prefix in MODEL_ORDER:
    f = find_file(prefix)
    if f is None:
        print(f"  [MISSING] {prefix:18s} -> file not found")
        datasets[prefix] = None
    else:
        print(f"  [OK]      {prefix:18s} -> {os.path.basename(f)}")
        datasets[prefix] = xr.open_dataset(f)

n_found = sum(1 for d in datasets.values() if d is not None)
print(f"\n  -> {n_found}/{len(MODEL_ORDER)} fichiers charges\n")

print("=" * 70)
print("GENERATION FIGURE MULTI-MODELE 5x2")
print("=" * 70)

fig = plt.figure(figsize=(18, 9))
gs = GridSpec(2, 5, figure=fig, hspace=0.28, wspace=0.10,
              left=0.03, right=0.98, top=0.91, bottom=0.16)

last_ims = None
for idx, prefix in enumerate(MODEL_ORDER):
    row, col = idx // 5, idx % 5
    ax = fig.add_subplot(gs[row, col], projection=ccrs.PlateCarree())
    label = MODEL_LABELS.get(prefix, prefix.upper())
    setup_subplot(ax, title=label)
    ds = datasets.get(prefix)
    if ds is not None:
        ims = plot_tercile(ax, ds)
        last_ims = ims
    else:
        ax.text(0.5, 0.5, f"File not found\n{prefix}",
                ha='center', va='center', transform=ax.transAxes,
                fontsize=10, color='red', style='italic')

if last_ims is not None:
    im_b, im_n, im_a = last_ims

    cax_b = fig.add_axes([0.08, 0.07, 0.24, 0.018])
    cb_b = plt.colorbar(im_b, cax=cax_b, orientation='horizontal', extend='max')
    cb_b.set_label('P(below-normal / dry) %', fontsize=10, fontweight='bold')
    cb_b.ax.tick_params(labelsize=9)
    cb_b.ax.invert_xaxis()

    cax_n = fig.add_axes([0.385, 0.07, 0.24, 0.018])
    cb_n = plt.colorbar(im_n, cax=cax_n, orientation='horizontal', extend='max')
    cb_n.set_label('P(near-normal) %', fontsize=10, fontweight='bold')
    cb_n.ax.tick_params(labelsize=9)

    cax_a = fig.add_axes([0.69, 0.07, 0.24, 0.018])
    cb_a = plt.colorbar(im_a, cax=cax_a, orientation='horizontal', extend='max')
    cb_a.set_label('P(above-normal / wet) %', fontsize=10, fontweight='bold')
    cb_a.ax.tick_params(labelsize=9)

fig.suptitle('Multi-System Precipitation Forecast — SON 2026 — North Africa\n'
              'Dominant tercile category per pixel (where probability >= 40%)',
              fontsize=15, fontweight='bold', y=0.97)

fig.text(0.5, 0.025,
          'C3S seasonal forecast systems (OSOP/SEAFOAM pipeline) — Initialization: August 2026  |  '
          'Calibration baseline: ERA5 1993-2016  |  '
          'Note: JMA CPS4 shown for reference but excluded from the MME.',
          ha='center', fontsize=8.5, style='italic', color='#555555')

out_file = os.path.join(OUT_DIR, "Fig4_OSOP_multisystem_TP_SON2026.png")
plt.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.2)
plt.close()
print(f"\n[OK] {out_file}")
print(f"     {os.path.getsize(out_file)/1024:.0f} KB")

for ds in datasets.values():
    if ds is not None:
        ds.close()

print("\n" + "=" * 70)
print(f"  explorer \"{OUT_DIR}\"")
