#!/usr/bin/env python
"""
Génération carte CCA NMME (Andrew Colman, UK Met Office) — v4
=============================================================
Correction v4 : un seul titre propre (suppression de la duplication)
"""
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import cptio

# ============================================================================
FC_DIR = r"C:/Users/wafae/Desktop/DATA PC HP/RCOFS/MEDCOF/MEDCOF26 JJA2026"
OUT_DIR = "C:/Users/wafae/projets/osop-data/seafoam/data/master/RCC_NA_maps"
os.makedirs(OUT_DIR, exist_ok=True)

EXTENT = [-18, 38, 14, 38]
THRESHOLD = 40.0
VAR = "precip"

MODELS = ["cola_ens","can4_ens","cfs_ens",
          "cola_sst","can4_sst","cfs_sst",
          "cola_ppn","can4_ppn","cfs_ppn"]
WEIGHTS = [2.0]*3 + [1.0]*6

# ============================================================================
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    return mcolors.LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n)))

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

CMAP_BELOW  = truncate_colormap(plt.cm.Blues, minval=0.3)
CMAP_NORMAL = truncate_colormap(plt.cm.Greys, minval=0.2, maxval=0.5)
CMAP_ABOVE  = plt.cm.YlOrRd

# ============================================================================
# 1. CHARGEMENT + MME
# ============================================================================
print("=" * 70)
print("MME NMME (Andrew Colman / UKMO)")
print("=" * 70)
weights = np.array(WEIGHTS); weights /= weights.sum()

mme = None
for model, w in zip(MODELS, weights):
    fcpath = os.path.join(FC_DIR, f"{model}.txt")
    print(f"  {model:12s}  poids={w:.4f}")
    ds = cptio.open_cptdataset(fcpath)
    if mme is None:
        mme = ds.copy(deep=True, data={VAR: np.zeros_like(ds[VAR].values)})
    mme[VAR].data += ds[VAR].data * w

print(f"\n✅ MME : {dict(mme.sizes)}")

# ============================================================================
# 2. MASQUAGE WMO
# ============================================================================
Z = mme[VAR].isel(T=0).values
Z1, Z2, Z3 = Z[0], Z[1], Z[2]

mask_below  = np.logical_or(np.logical_or(Z1 < Z2, Z1 < Z3), Z1 < THRESHOLD)
mask_normal = np.logical_or(np.logical_or(Z2 < Z1, Z2 < Z3), Z2 < THRESHOLD)
mask_above  = np.logical_or(np.logical_or(Z3 < Z1, Z3 < Z2), Z3 < THRESHOLD)

Z_below_show  = np.where(~mask_below,  Z1, np.nan)
Z_normal_show = np.where(~mask_normal, Z2, np.nan)
Z_above_show  = np.where(~mask_above,  Z3, np.nan)

lat = mme.coords['Y'].values
lon = mme.coords['X'].values

# ============================================================================
# 3. CARTE — UN SEUL titre propre
# ============================================================================
print("\n" + "=" * 70)
print("Génération Fig 13 v4...")
print("=" * 70)

fig = plt.figure(figsize=(12, 7))
ax = fig.add_axes([0.05, 0.22, 0.90, 0.62], projection=ccrs.PlateCarree())
ax.set_extent(EXTENT, crs=ccrs.PlateCarree())

# Données (zorder bas)
levels = np.arange(THRESHOLD, 85.0, 5.0)
norm_below  = BoundaryNorm(levels, ncolors=CMAP_BELOW.N, extend="max")
norm_normal = BoundaryNorm(levels, ncolors=CMAP_NORMAL.N)
norm_above  = BoundaryNorm(levels, ncolors=CMAP_ABOVE.N, extend="max")

im_below = ax.pcolormesh(lon, lat, Z_below_show, transform=ccrs.PlateCarree(),
                          cmap=CMAP_BELOW, norm=norm_below, shading='auto', zorder=2)
im_normal = ax.pcolormesh(lon, lat, Z_normal_show, transform=ccrs.PlateCarree(),
                           cmap=CMAP_NORMAL, norm=norm_normal, shading='auto', zorder=2)
im_above = ax.pcolormesh(lon, lat, Z_above_show, transform=ccrs.PlateCarree(),
                          cmap=CMAP_ABOVE, norm=norm_above, shading='auto', zorder=2)

# Océan PAR-DESSUS les données pour masquer le débordement
ax.add_feature(cfeature.OCEAN, facecolor='#E8F4F8', zorder=10)
ax.add_feature(cfeature.LAKES, facecolor='#E8F4F8', zorder=10)
ax.add_feature(cfeature.COASTLINE, linewidth=0.7, edgecolor='black', zorder=11)
ax.add_feature(MA_BORDERS, linewidth=0.5, edgecolor='black', zorder=11)

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray',
                   linestyle='--', alpha=0.5, zorder=12)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8}
gl.ylabel_style = {'size': 8}

# ✨ UN SEUL TITRE — clean et complet (au lieu de suptitle + set_title)
# Pas de ax.set_title() ni fig.suptitle() qui se dupliquent
fig.text(0.5, 0.92,
          'Figure 13 — NMME CCA Multi-Method T2m Forecast for JJA 2026',
          ha='center', fontsize=13, fontweight='bold')
fig.text(0.5, 0.88,
          'Dominant tercile category per pixel (where probability ≥ 40%)',
          ha='center', fontsize=10, style='italic', color='#444444')

# 3 colorbars
cax1 = fig.add_axes([0.08, 0.10, 0.26, 0.028])
cb1 = plt.colorbar(im_below, cax=cax1, orientation='horizontal', extend='max')
cb1.set_label('P(below) %', fontsize=9, fontweight='bold', labelpad=4)
cb1.ax.tick_params(labelsize=8)
cb1.ax.invert_xaxis()

cax2 = fig.add_axes([0.385, 0.10, 0.26, 0.028])
cb2 = plt.colorbar(im_normal, cax=cax2, orientation='horizontal', extend='neither')
cb2.set_label('P(normal) %', fontsize=9, fontweight='bold', labelpad=4)
cb2.ax.tick_params(labelsize=8)

cax3 = fig.add_axes([0.69, 0.10, 0.26, 0.028])
cb3 = plt.colorbar(im_above, cax=cax3, orientation='horizontal', extend='max')
cb3.set_label('P(above) %', fontsize=9, fontweight='bold', labelpad=4)
cb3.ax.tick_params(labelsize=8)

# Note méthodologique tout en bas
fig.text(0.5, 0.02,
         '3 GCMs (COLA, CFS, CanSIPSv3) × 3 methods: Direct ensemble (×2) + CCA SST (×1) + CCA precip (×1)  |  '
         'Initialization: May 2026  |  Source: UK Met Office',
         ha='center', fontsize=8, style='italic', color='#555555')

out_file = os.path.join(OUT_DIR, "Fig13_NMME_CCA_T2m_JJA2026.png")
plt.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.2)
plt.close()
print(f"\n✅ {out_file}")
print(f"   {os.path.getsize(out_file)/1024:.0f} KB")
print("\n" + "=" * 70)
print(f"  explorer \"{OUT_DIR}\"")
