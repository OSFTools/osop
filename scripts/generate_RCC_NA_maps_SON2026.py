#!/usr/bin/env python
"""
Génération des cartes RCC North Africa — SON 2026 — Version v5
================================================================
STYLE OSOP/WMO : une seule carte par variable montre les 3 terciles
combinés (below / normal / above), où chaque pixel est colorié selon
le tercile dominant. 3 colorbars côte à côte.

Cartes produites :
  Fig 5 : Skill Matrix (4 panels Pearson + ROCSS) — palette OSOP RdYlBu_r
  Fig 6 : Forecast unique T2m (3 terciles superposés)
  Fig 7 : Forecast unique TP (3 terciles superposés, VERT pour humide)
  Fig 8 : Sub-regional synthesis (avec annotations pays)
"""
import os
import glob
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE = 'C:/Users/wafae/projets/osop-data/seafoam/data/master'
HIND_SCORES = f'{BASE}/hindcast/scores'
FCST_PRODUCTS = f'{BASE}/forecast/products'
OUT_DIR = f'{BASE}/RCC_NA_maps'
os.makedirs(OUT_DIR, exist_ok=True)

EXTENT = [-18, 38, 14, 38]

# Seuil WMO : tercile dominant exploité seulement si > 40%
THRESHOLD = 40.0

# ============================================================================
# FRONTIÈRES MAROC OFFICIELLES (méthode OSOP)
# ============================================================================
def get_morocco_borders():
    border_set = "admin_0_countries_mar"
    try:
        _ = shpreader.natural_earth(resolution="10m", category="cultural", name=border_set)
    except KeyError:
        _ = shpreader.natural_earth(resolution="10m", category="cultural", name=border_set)
    finally:
        local = cfeature.NaturalEarthFeature(
            category="cultural", name=border_set, scale="10m", facecolor="none"
        )
    return local

print("Chargement du shapefile Maroc (admin_0_countries_mar)...")
MA_BORDERS = get_morocco_borders()
print("[OK] Shapefile chargé\n")

# ============================================================================
# UTILS
# ============================================================================
def find(pattern, directory):
    f = glob.glob(os.path.join(directory, pattern))
    return f[0] if f else None

def setup_ax(ax, title='', extent=EXTENT):
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6, edgecolor='black')
    ax.add_feature(MA_BORDERS, linewidth=0.5, edgecolor='black')
    ax.add_feature(cfeature.OCEAN, facecolor='#E8F4F8', zorder=0)
    ax.add_feature(cfeature.LAND,  facecolor='#FAFAF5', zorder=0)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 7}
    gl.ylabel_style = {'size': 7}
    if title:
        ax.set_title(title, fontsize=10, fontweight='bold', pad=8)

# ============================================================================
# FONCTION CLÉ : carte tercile unique (style OSOP/WMO)
# ============================================================================
def plot_tercile_map(ax, lon, lat, lower, middle, higher,
                      cmap_below, cmap_above, cmap_normal,
                      threshold=THRESHOLD,
                      levels_high=None, levels_normal=None):
    """
    Affiche une carte montrant le tercile dominant en couleur,
    et renvoie les 3 objets pcolormesh (pour les 3 colorbars).

    Méthode :
      1. Pour chaque pixel, déterminer la catégorie gagnante (max des 3 terciles)
      2. Masquer les zones où la proba dominante < threshold (=> apparaît en blanc/gris)
      3. Tracer 3 couches superposées avec masques disjoints
    """
    if levels_high is None:
        levels_high = np.arange(40, 95, 5)  # 40..90 par 5
    if levels_normal is None:
        levels_normal = np.arange(40, 55, 5)  # 40, 45, 50

    # Déterminer le tercile dominant par pixel
    stack = np.stack([lower, middle, higher], axis=0)  # (3, lat, lon)
    dominant = np.argmax(stack, axis=0)  # 0=below, 1=normal, 2=above

    # Masques disjoints (tercile dominant ET au-dessus du seuil)
    mask_below  = (dominant == 0) & (lower  >= threshold)
    mask_normal = (dominant == 1) & (middle >= threshold)
    mask_above  = (dominant == 2) & (higher >= threshold)

    # Fields masqués pour affichage
    field_below  = np.where(mask_below,  lower,  np.nan)
    field_normal = np.where(mask_normal, middle, np.nan)
    field_above  = np.where(mask_above,  higher, np.nan)

    # Tracé des 3 couches
    norm_high = BoundaryNorm(np.arange(40, 95, 5), ncolors=256, extend='max')
    norm_normal = BoundaryNorm(np.arange(40, 55, 5), ncolors=256, extend='max')

    im_below = ax.pcolormesh(lon, lat, field_below,
                              transform=ccrs.PlateCarree(),
                              cmap=cmap_below, norm=norm_high, shading='auto')
    im_normal = ax.pcolormesh(lon, lat, field_normal,
                               transform=ccrs.PlateCarree(),
                               cmap=cmap_normal, norm=norm_normal, shading='auto')
    im_above = ax.pcolormesh(lon, lat, field_above,
                              transform=ccrs.PlateCarree(),
                              cmap=cmap_above, norm=norm_high, shading='auto')
    return im_below, im_normal, im_above

def add_3_colorbars(fig, im_below, im_normal, im_above,
                     label_below='Prob. below (%)',
                     label_normal='Prob. normal (%)',
                     label_above='Prob. above (%)'):
    """Ajoute 3 colorbars horizontales côte à côte sous la carte."""
    cax1 = fig.add_axes([0.13, 0.05, 0.22, 0.025])
    cb1 = plt.colorbar(im_below, cax=cax1, orientation='horizontal', extend='max')
    cb1.set_label(label_below, fontsize=9, fontweight='bold')
    cb1.ax.tick_params(labelsize=8)
    cb1.ax.invert_xaxis()  # convention OSOP : below = inversé

    cax2 = fig.add_axes([0.4, 0.05, 0.22, 0.025])
    cb2 = plt.colorbar(im_normal, cax=cax2, orientation='horizontal', extend='max')
    cb2.set_label(label_normal, fontsize=9, fontweight='bold')
    cb2.ax.tick_params(labelsize=8)

    cax3 = fig.add_axes([0.67, 0.05, 0.22, 0.025])
    cb3 = plt.colorbar(im_above, cax=cax3, orientation='horizontal', extend='max')
    cb3.set_label(label_above, fontsize=9, fontweight='bold')
    cb3.ax.tick_params(labelsize=8)

# ============================================================================
# PALETTES OSOP/WMO
# ============================================================================
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    return mcolors.LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n))
    )

# Pour T2m (Bleu / Gris / Rouge-orange)
T2M_CMAP_BELOW  = truncate_colormap(plt.cm.Blues,  minval=0.3)
T2M_CMAP_NORMAL = truncate_colormap(plt.cm.Greys,  minval=0.2, maxval=0.5)
T2M_CMAP_ABOVE  = plt.cm.YlOrRd

# Pour TP (Brun / Gris / VERT)
TP_CMAP_BELOW  = truncate_colormap(plt.cm.YlOrBr,  minval=0.3)  # Brun
TP_CMAP_NORMAL = truncate_colormap(plt.cm.Greys,   minval=0.2, maxval=0.5)
TP_CMAP_ABOVE  = truncate_colormap(plt.cm.Greens,  minval=0.3)  # VERT

# ============================================================================
# CARTE 5 — SKILL MATRIX (4 panels — palette OSOP RdYlBu_r)
# ============================================================================
print("=" * 70)
print("CARTE 5 — Skill Matrix")
print("=" * 70)

fig = plt.figure(figsize=(14, 9))
gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.15)

score_levels = np.linspace(-1.0, 1.0, 11)

# (a) Pearson T2m
ax1 = fig.add_subplot(gs[0,0], projection=ccrs.PlateCarree())
setup_ax(ax1, '(a) MME Pearson correlation — T2m SON (hindcast 1993–2016)')
f = find('mme_1_1993-2016*8_234*2m_temperature*3m*pearson_corr.nc', HIND_SCORES)
ds = xr.open_dataset(f); var = list(ds.data_vars)[0]
im = ax1.contourf(ds.lon, ds.lat, ds[var].isel(forecastMonth=0),
                   levels=score_levels, transform=ccrs.PlateCarree(),
                   cmap='RdYlBu_r', extend='both')
plt.colorbar(im, ax=ax1, orientation='horizontal', shrink=0.85, pad=0.06,
              label='Pearson r')
ds.close()

# (b) Pearson TP
ax2 = fig.add_subplot(gs[0,1], projection=ccrs.PlateCarree())
setup_ax(ax2, '(b) MME Pearson correlation — TP SON (hindcast 1993–2016)')
f = find('mme_1_1993-2016*8_234*total_precipitation*3m*pearson_corr.nc', HIND_SCORES)
ds = xr.open_dataset(f); var = list(ds.data_vars)[0]
im = ax2.contourf(ds.lon, ds.lat, ds[var].isel(forecastMonth=0),
                   levels=score_levels, transform=ccrs.PlateCarree(),
                   cmap='RdYlBu_r', extend='both')
plt.colorbar(im, ax=ax2, orientation='horizontal', shrink=0.85, pad=0.06,
              label='Pearson r')
ds.close()

# (c) ROCSS T2m above
ax3 = fig.add_subplot(gs[1,0], projection=ccrs.PlateCarree())
setup_ax(ax3, '(c) MME ROCSS — T2m above-normal SON')
f = find('mme_1_1993-2016*8_234*2m_temperature*3m*rocss.nc', HIND_SCORES)
ds = xr.open_dataset(f); var = list(ds.data_vars)[0]
im = ax3.contourf(ds.lon, ds.lat,
                   ds[var].isel(forecastMonth=0, category=2),
                   levels=score_levels, transform=ccrs.PlateCarree(),
                   cmap='RdYlBu_r', extend='both')
plt.colorbar(im, ax=ax3, orientation='horizontal', shrink=0.85, pad=0.06,
              label='ROCSS')
ds.close()

# (d) ROCSS TP below
ax4 = fig.add_subplot(gs[1,1], projection=ccrs.PlateCarree())
setup_ax(ax4, '(d) MME ROCSS — TP below-normal SON')
f = find('mme_1_1993-2016*8_234*total_precipitation*3m*rocss.nc', HIND_SCORES)
ds = xr.open_dataset(f); var = list(ds.data_vars)[0]
im = ax4.contourf(ds.lon, ds.lat,
                   ds[var].isel(forecastMonth=0, category=0),
                   levels=score_levels, transform=ccrs.PlateCarree(),
                   cmap='RdYlBu_r', extend='both')
plt.colorbar(im, ax=ax4, orientation='horizontal', shrink=0.85, pad=0.06,
              label='ROCSS')
ds.close()

fig.suptitle('Figure 5 — MME Forecast Skill over North Africa\nHindcast verification 1993–2016 vs ERA5',
             fontsize=13, fontweight='bold', y=0.995)
out5 = f'{OUT_DIR}/Fig5_MME_skill_matrix.png'
plt.savefig(out5, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] {out5}")

# ============================================================================
# CARTE 6 — FORECAST T2m (3 terciles superposés)
# ============================================================================
print("\n" + "=" * 70)
print("CARTE 6 — Forecast T2m (3 tercile categories combined)")
print("=" * 70)

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
setup_ax(ax,
          'Probabilistic Forecast — T2m SON 2026\nDominant tercile category per pixel (where probability ≥ 40%)')

f = find('mme_1_2026-2026*8_234*2m_temperature*3m*forecast_percentages.nc', FCST_PRODUCTS)
ds = xr.open_dataset(f)
lower  = ds['lower'].isel(forecastMonth=0).values
middle = ds['middle'].isel(forecastMonth=0).values
higher = ds['higher'].isel(forecastMonth=0).values

im_below, im_normal, im_above = plot_tercile_map(
    ax, ds.lon, ds.lat, lower, middle, higher,
    cmap_below=T2M_CMAP_BELOW,
    cmap_above=T2M_CMAP_ABOVE,
    cmap_normal=T2M_CMAP_NORMAL
)
ds.close()

# Ajustement layout pour laisser place aux 3 colorbars
plt.subplots_adjust(bottom=0.15, top=0.92)
add_3_colorbars(fig, im_below, im_normal, im_above,
                 label_below='P(below) %',
                 label_normal='P(normal) %',
                 label_above='P(above) %')

fig.suptitle('Figure 6 — MME Forecast for T2m — SON 2026',
             fontsize=13, fontweight='bold', y=0.99)
out6 = f'{OUT_DIR}/Fig6_MME_forecast_T2m_SON2026.png'
plt.savefig(out6, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] {out6}")

# ============================================================================
# CARTE 7 — FORECAST TP (3 terciles superposés, VERT pour humide)
# ============================================================================
print("\n" + "=" * 70)
print("CARTE 7 — Forecast TP (3 tercile categories, GREEN for wet)")
print("=" * 70)

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
setup_ax(ax,
          'Probabilistic Forecast — TP SON 2026\nDominant tercile category per pixel (where probability ≥ 40%)')

f = find('mme_1_2026-2026*8_234*total_precipitation*3m*forecast_percentages.nc', FCST_PRODUCTS)
ds = xr.open_dataset(f)
lower  = ds['lower'].isel(forecastMonth=0).values
middle = ds['middle'].isel(forecastMonth=0).values
higher = ds['higher'].isel(forecastMonth=0).values

im_below, im_normal, im_above = plot_tercile_map(
    ax, ds.lon, ds.lat, lower, middle, higher,
    cmap_below=TP_CMAP_BELOW,   # brun = sec
    cmap_above=TP_CMAP_ABOVE,   # vert = humide
    cmap_normal=TP_CMAP_NORMAL
)
ds.close()

plt.subplots_adjust(bottom=0.15, top=0.92)
add_3_colorbars(fig, im_below, im_normal, im_above,
                 label_below='P(below/dry) %',
                 label_normal='P(normal) %',
                 label_above='P(above/wet) %')

fig.suptitle('Figure 7 — MME Forecast for Precipitation — SON 2026',
             fontsize=13, fontweight='bold', y=0.99)
out7 = f'{OUT_DIR}/Fig7_MME_forecast_TP_SON2026.png'
plt.savefig(out7, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] {out7}")

# ============================================================================
# CARTE 8 — SUB-REGIONAL SYNTHESIS (2 panels avec annotations)
# ============================================================================
print("\n" + "=" * 70)
print("CARTE 8 — Sub-Regional Synthesis (with country annotations)")
print("=" * 70)

COUNTRY_BOXES = {
    'Morocco':    {'box': (-13.0, 27.5, -1.0, 35.9), 't2m_above': 75.4, 'tp_below': 42.2},
    'Algeria N':  {'box': ( -8.0, 32.0, 12.0, 37.0), 't2m_above': 79.4, 'tp_below': 46.7},
    'Algeria S':  {'box': ( -5.0, 19.0, 11.0, 28.0), 't2m_above': 82.2, 'tp_below': 35.7},
    'Tunisia':    {'box': (  7.5, 30.0, 11.6, 37.4), 't2m_above': 79.7, 'tp_below': 45.1},
    'Libya':      {'box': (  9.4, 19.5, 25.2, 33.2), 't2m_above': 77.3, 'tp_below': 33.7},
    'Egypt':      {'box': ( 24.7, 22.0, 36.9, 31.7), 't2m_above': 85.0, 'tp_below': 23.6},
    'Mauritania': {'box': (-17.1, 14.7, -4.8, 27.3), 't2m_above': 57.0, 'tp_below': 44.9},
}

fig = plt.figure(figsize=(14, 8))
gs = GridSpec(1, 2, figure=fig, wspace=0.08)

# (a) T2m signal
ax1 = fig.add_subplot(gs[0,0], projection=ccrs.PlateCarree())
setup_ax(ax1, '(a) T2m Signal — SON 2026 (MME, dominant tercile)')
f = find('mme_1_2026-2026*8_234*2m_temperature*3m*forecast_percentages.nc', FCST_PRODUCTS)
ds = xr.open_dataset(f)
lower  = ds['lower'].isel(forecastMonth=0).values
middle = ds['middle'].isel(forecastMonth=0).values
higher = ds['higher'].isel(forecastMonth=0).values
im_b1, im_n1, im_a1 = plot_tercile_map(
    ax1, ds.lon, ds.lat, lower, middle, higher,
    cmap_below=T2M_CMAP_BELOW, cmap_above=T2M_CMAP_ABOVE,
    cmap_normal=T2M_CMAP_NORMAL)
ds.close()
for country, info in COUNTRY_BOXES.items():
    lon_min, lat_min, lon_max, lat_max = info['box']
    rect = Rectangle((lon_min, lat_min), lon_max-lon_min, lat_max-lat_min,
                      linewidth=1.5, edgecolor='black', facecolor='none',
                      transform=ccrs.PlateCarree())
    ax1.add_patch(rect)
    cx, cy = (lon_min+lon_max)/2, (lat_min+lat_max)/2
    ax1.text(cx, cy, f"{country}\nT2m+: {info['t2m_above']}%",
              fontsize=8, fontweight='bold', ha='center', va='center',
              transform=ccrs.PlateCarree(),
              bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       edgecolor='black', alpha=0.85))

# (b) TP signal
ax2 = fig.add_subplot(gs[0,1], projection=ccrs.PlateCarree())
setup_ax(ax2, '(b) TP Signal — SON 2026 (MME, dominant tercile)')
f = find('mme_1_2026-2026*8_234*total_precipitation*3m*forecast_percentages.nc', FCST_PRODUCTS)
ds = xr.open_dataset(f)
lower  = ds['lower'].isel(forecastMonth=0).values
middle = ds['middle'].isel(forecastMonth=0).values
higher = ds['higher'].isel(forecastMonth=0).values
im_b2, im_n2, im_a2 = plot_tercile_map(
    ax2, ds.lon, ds.lat, lower, middle, higher,
    cmap_below=TP_CMAP_BELOW, cmap_above=TP_CMAP_ABOVE,
    cmap_normal=TP_CMAP_NORMAL)
ds.close()
for country, info in COUNTRY_BOXES.items():
    lon_min, lat_min, lon_max, lat_max = info['box']
    rect = Rectangle((lon_min, lat_min), lon_max-lon_min, lat_max-lat_min,
                      linewidth=1.5, edgecolor='black', facecolor='none',
                      transform=ccrs.PlateCarree())
    ax2.add_patch(rect)
    cx, cy = (lon_min+lon_max)/2, (lat_min+lat_max)/2
    ax2.text(cx, cy, f"{country}\nTP-: {info['tp_below']}%",
              fontsize=8, fontweight='bold', ha='center', va='center',
              transform=ccrs.PlateCarree(),
              bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       edgecolor='black', alpha=0.85))

fig.suptitle('Figure 8 — Sub-Regional Synthesis SON 2026\nMME tercile probabilities per country/region',
             fontsize=13, fontweight='bold', y=1.02)
out8 = f'{OUT_DIR}/Fig8_summary_subregions.png'
plt.savefig(out8, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] {out8}")

# ============================================================================
print("\n" + "=" * 70)
print("[OK] CARTES RÉGÉNÉRÉES (style OSOP — tercile dominant)")
print("=" * 70)
print(f"\nDossier : {OUT_DIR}")
for f in sorted(os.listdir(OUT_DIR)):
    full = os.path.join(OUT_DIR, f)
    print(f"  • {f}  ({os.path.getsize(full)/1024:.0f} KB)")
print(f"\n  explorer \"{OUT_DIR}\"\n")
