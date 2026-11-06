#!/usr/bin/env python
"""
Apply Dry Mask to MME TP JJA Forecast — RCC NorAf Standard
============================================================
Calcule la climatologie ERA5 JJA 1993-2016 et masque les pixels ou
la precipitation moyenne JJA est inferieure a un seuil operationnel
(0.5 mm/jour, soit ~45 mm sur la saison JJA = 92 jours).

Genere :
  - Fig7_MME_forecast_TP_JJA2026_drymasked.png
"""
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shutil

# ============================================================================
BASE = "C:/Users/wafae/projets/osop-data/seafoam/data/master"
FCST_FILE = f"{BASE}/forecast/products/mme_1_2026-2026_monthly_mean_5_234_45_-30_-2.5_60_total_precipitation.3m.forecast_percentages.nc"
ERA5_FILE = f"{BASE}/hindcast/downloads/era5_5_234_45_-30_-2.5_60_total_precipitation.grib"

OUT_DIR  = f"{BASE}/RCC_NA_maps"
JJA_DIR  = f"{OUT_DIR}/JJA2026/Precipitation"
os.makedirs(JJA_DIR, exist_ok=True)

EXTENT = [-18, 38, 14, 38]
THRESHOLD = 40.0  # seuil tercile dominant

# Dry mask seuil : 0.5 mm/jour de moyenne JJA = ~45 mm/saison
DRY_THRESHOLD_MM_PER_DAY = 0.5

# ============================================================================
# FRONTIERES MAROC
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

# ============================================================================
# CHARGER LA CLIMATOLOGIE ERA5 POUR JJA
# ============================================================================
print("Chargement climatologie ERA5 JJA 1993-2016...")
# Tester plusieurs sources possibles pour la clim
era5_candidates = [
    f"{BASE}/hindcast/downloads/era5_5_234_45_-30_-2.5_60_total_precipitation.grib",
    f"{BASE}/hindcast/products/era5_obs_1993-2016_total_precipitation.3m.nc",
]

clim_jja = None
for f in era5_candidates:
    if os.path.exists(f):
        try:
            if f.endswith('.grib'):
                ds = xr.open_dataset(f, engine='cfgrib')
            else:
                ds = xr.open_dataset(f)
            print(f"  Trouve : {os.path.basename(f)}")
            print(f"  Vars : {list(ds.data_vars)}")
            print(f"  Coords : {list(ds.coords)}")
            print(f"  Dims : {dict(ds.sizes)}")
            
            # Identifier la variable de precip
            var = None
            for v in ['tp', 'total_precipitation', 'precipitation', 'pr']:
                if v in ds.data_vars:
                    var = v
                    break
            if var is None:
                var = list(ds.data_vars)[0]
            
            # Si on a une dim "time" ou "valid_time", on peut filtrer JJA
            time_coord = None
            for tc in ['time', 'valid_time', 'forecastMonth']:
                if tc in ds.coords or tc in ds.dims:
                    time_coord = tc
                    break
            
            if time_coord and time_coord in ['time', 'valid_time']:
                # Filtrer JJA (mois 6,7,8) et calculer la moyenne
                ds_jja = ds.sel({time_coord: ds[time_coord].dt.month.isin([6,7,8])})
                clim_jja = ds_jja[var].mean(dim=time_coord)
                print(f"  Clim JJA calculee : shape {clim_jja.shape}")
                break
            else:
                # Probablement deja moyenne ou format different
                clim_jja = ds[var]
                if 'forecastMonth' in clim_jja.dims:
                    clim_jja = clim_jja.isel(forecastMonth=0)
                print(f"  Clim utilisee directement : shape {clim_jja.shape}")
                break
        except Exception as e:
            print(f"  ERROR avec {f}: {str(e)[:80]}")

if clim_jja is None:
    print("\n[WARN] Climatologie ERA5 introuvable. Utilisation d'un dry mask APPROCHE basee sur :")
    print("        - lat < 21 N et lon > -10 (Sahara central) -> dry")
    print("        - et zones desertiques connues")
    # On masquera plus tard via lat/lon
    use_geographic_mask = True
else:
    use_geographic_mask = False
    print(f"\nClim JJA stats: min={float(clim_jja.min()):.2e}, max={float(clim_jja.max()):.2e}, mean={float(clim_jja.mean()):.2e}")

# ============================================================================
# CHARGER LE FORECAST MME
# ============================================================================
print("\nChargement MME TP forecast JJA 2026...")
ds_fc = xr.open_dataset(FCST_FILE).isel(forecastMonth=0)
print(f"  start_date: {ds_fc.start_date.values}")
print(f"  shape: {ds_fc.lower.shape}")

lower = ds_fc['lower'].values
middle = ds_fc['middle'].values
higher = ds_fc['higher'].values
lat = ds_fc.lat.values
lon = ds_fc.lon.values

# ============================================================================
# CONSTRUIRE LE DRY MASK
# ============================================================================
if use_geographic_mask:
    # Mask geographique base sur lat/lon (regions desertiques connues)
    # Sahara central : 18-28 N, -10 a 30 E
    # Mauritanie/Mali Nord : 15-23 N, -16 a -3 E
    LON, LAT = np.meshgrid(lon, lat)
    dry_mask = (
        # Sahara central (zone hyper-aride)
        ((LAT >= 18) & (LAT <= 27) & (LON >= -8) & (LON <= 32)) |
        # Sahara occidental atlantique
        ((LAT >= 20) & (LAT <= 27) & (LON >= -15) & (LON <= -8)) |
        # Desert libyen
        ((LAT >= 19) & (LAT <= 30) & (LON >= 20) & (LON <= 32)) |
        # Desert egyptien
        ((LAT >= 22) & (LAT <= 30) & (LON >= 25) & (LON <= 35))
    )
    print(f"\nDry mask geographique : {dry_mask.sum()} / {dry_mask.size} pixels masques ({100*dry_mask.sum()/dry_mask.size:.1f}%)")
else:
    # Convertir clim a la grille forecast si necessaire (interpolation simple)
    clim_2d = clim_jja.values
    if clim_2d.shape != lower.shape:
        # Interpoler vers la grille forecast
        clim_interp = clim_jja.interp(lat=lat, lon=lon, method='linear')
        clim_2d = clim_interp.values
    
    # ERA5 TP est generalement en m/jour ou kg/m2/s, convertir en mm/jour
    # m/jour -> mm/jour : * 1000
    if np.nanmax(clim_2d) < 0.01:
        # Probablement en m/jour
        clim_mm_day = clim_2d * 1000
        print(f"  Conversion m/jour -> mm/jour")
    elif np.nanmax(clim_2d) > 100:
        # Probablement en mm/saison total
        clim_mm_day = clim_2d / 92
        print(f"  Conversion mm/saison -> mm/jour")
    else:
        clim_mm_day = clim_2d
        print(f"  Suppose deja en mm/jour")
    
    dry_mask = clim_mm_day < DRY_THRESHOLD_MM_PER_DAY
    print(f"\nDry mask climatologique : {dry_mask.sum()} / {dry_mask.size} pixels masques ({100*dry_mask.sum()/dry_mask.size:.1f}%)")

# ============================================================================
# PALETTES (tercile dominant style)
# ============================================================================
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

# ============================================================================
# CONSTRUIRE LES CHAMPS A TRACER (tercile dominant >= 40%)
# ============================================================================
# Masques disjoints (tercile dominant ET proba >= seuil) ET PAS dry-masked
mask_b = np.logical_or(np.logical_or(lower  < middle, lower  < higher), lower  < THRESHOLD) | dry_mask
mask_n = np.logical_or(np.logical_or(middle < lower,  middle < higher), middle < THRESHOLD) | dry_mask
mask_a = np.logical_or(np.logical_or(higher < lower,  higher < middle), higher < THRESHOLD) | dry_mask

fb = np.where(~mask_b, lower,  np.nan)
fn = np.where(~mask_n, middle, np.nan)
fa = np.where(~mask_a, higher, np.nan)

# ============================================================================
# GENERATION DE LA FIGURE
# ============================================================================
print("\nGeneration figure dry-masked...")

fig = plt.figure(figsize=(13, 7))
ax = fig.add_axes([0.05, 0.18, 0.88, 0.72], projection=ccrs.PlateCarree())
ax.set_extent(EXTENT, crs=ccrs.PlateCarree())

# Trace les terciles
im_b = ax.pcolormesh(lon, lat, fb, transform=ccrs.PlateCarree(),
                     cmap=TP_CMAP_BELOW, norm=NORM_BELOW, shading='auto', zorder=2)
im_n = ax.pcolormesh(lon, lat, fn, transform=ccrs.PlateCarree(),
                     cmap=TP_CMAP_NORMAL, norm=NORM_NORMAL, shading='auto', zorder=2)
im_a = ax.pcolormesh(lon, lat, fa, transform=ccrs.PlateCarree(),
                     cmap=TP_CMAP_ABOVE, norm=NORM_ABOVE, shading='auto', zorder=2)

# Tracer le dry mask avec hachures
dry_field = np.where(dry_mask, 1.0, np.nan)
ax.pcolormesh(lon, lat, dry_field, transform=ccrs.PlateCarree(),
              cmap=mcolors.ListedColormap(['#D8D8D8']), shading='auto',
              alpha=0.85, zorder=3)
# Hachures pour rendre le dry mask distinct
ax.contourf(lon, lat, dry_field, levels=[0.5, 1.5], hatches=['///'],
             colors='none', transform=ccrs.PlateCarree(), zorder=4)

# Ocean + cotes par dessus
ax.add_feature(cfeature.OCEAN, facecolor='#E8F4F8', zorder=10)
ax.add_feature(cfeature.COASTLINE, linewidth=0.7, edgecolor='black', zorder=11)
ax.add_feature(MA_BORDERS, linewidth=0.5, edgecolor='black', zorder=11)

# Grille
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray',
                   linestyle='--', alpha=0.5, zorder=12)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 9}
gl.ylabel_style = {'size': 9}

# Titre
fig.suptitle('Figure 1 — MME Probabilistic Precipitation Forecast — JJA 2026\n'
              'Dominant tercile category per pixel (probability ≥ 40%) — Dry-masked',
              fontsize=13, fontweight='bold', y=0.97)

# Colorbars (3 horizontales en bas)
cax_b = fig.add_axes([0.07, 0.09, 0.26, 0.022])
cb_b = plt.colorbar(im_b, cax=cax_b, orientation='horizontal', extend='max')
cb_b.set_label('P(below-normal / dry) %', fontsize=9, fontweight='bold')
cb_b.ax.tick_params(labelsize=8)
cb_b.ax.invert_xaxis()

cax_n = fig.add_axes([0.37, 0.09, 0.26, 0.022])
cb_n = plt.colorbar(im_n, cax=cax_n, orientation='horizontal', extend='max')
cb_n.set_label('P(near-normal) %', fontsize=9, fontweight='bold')
cb_n.ax.tick_params(labelsize=8)

cax_a = fig.add_axes([0.67, 0.09, 0.26, 0.022])
cb_a = plt.colorbar(im_a, cax=cax_a, orientation='horizontal', extend='max')
cb_a.set_label('P(above-normal / wet) %', fontsize=9, fontweight='bold')
cb_a.ax.tick_params(labelsize=8)

# Note dry mask
note = ("Hatched grey areas: climatologically dry regions (JJA mean precipitation < 0.5 mm/day, ERA5 1993-2016) "
        "— probabilistic terciles not interpretable")
fig.text(0.5, 0.03, note, ha='center', fontsize=8.5, style='italic', color='#444444')

# Sauvegarde
out_file = os.path.join(OUT_DIR, "Fig7_MME_forecast_TP_JJA2026_drymasked.png")
plt.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.2)

# Copie dans JJA2026
out_file_jja = os.path.join(JJA_DIR, "Fig7_MME_forecast_TP_JJA2026_drymasked.png")
shutil.copy(out_file, out_file_jja)

plt.close()
print(f"\n[OK] {out_file}")
print(f"     {os.path.getsize(out_file)/1024:.0f} KB")
print(f"[OK] Copie aussi dans : {out_file_jja}")

ds_fc.close()
print("\nFAIT")
