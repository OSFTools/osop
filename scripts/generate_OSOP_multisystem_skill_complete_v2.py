#!/usr/bin/env python
"""
OSOP/SEAFOAM Multi-System Skill Matrix v2 — JJA hindcast 1993-2016 — North Africa
==================================================================================
Genere pour CHAQUE score restant et CHAQUE variable (TP, T2m) une figure 10 panels :
  - spearman_corr            -> 1 figure
  - bs (Brier Score above)   -> 1 figure
  - rps (Ranked Prob Score)  -> 1 figure
  - roc (ROC area above)     -> 1 figure
  - rocss below              -> 1 figure
Total : 5 figures x 2 variables = 10 figures

Toutes copiees automatiquement dans JJA2026/Precipitation/ et JJA2026/TMP/

Convention category : [0=below, 1=normal, 2=above]
"""
import os
import glob
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shutil

# ============================================================================
SCORES_DIR = "C:/Users/wafae/projets/osop-data/seafoam/data/master/hindcast/scores"
OUT_DIR    = "C:/Users/wafae/projets/osop-data/seafoam/data/master/RCC_NA_maps"
JJA_DIR    = os.path.join(OUT_DIR, "JJA2026")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(os.path.join(JJA_DIR, "TMP"), exist_ok=True)
os.makedirs(os.path.join(JJA_DIR, "Precipitation"), exist_ok=True)

EXTENT = [-18, 38, 14, 38]

MODEL_LABELS = {
    'ecmwf_51':       'ECMWF SEAS5',
    'ukmo_610':       'UK Met Office GS6.1',
    'meteo_france_9': 'Meteo-France Sys 9',
    'dwd_22':         'DWD GCFS 2.1',
    'cmcc_4':         'CMCC SPS 4.0',
    'ncep_2':         'NCEP CFSv2',
    'eccc_4':         'ECCC CanSIPSv3',
    'eccc_5':         'ECCC GEM5',
    'jma_4':          'JMA CPS4',
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

# ============================================================================
def find_score_file(model_prefix, variable, score_name):
    pattern = f"{model_prefix}_1993-2016*{variable}*3m.{score_name}.nc"
    matches = glob.glob(os.path.join(SCORES_DIR, pattern))
    return matches[0] if matches else None

def setup_subplot(ax, title=''):
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, facecolor='#E8F4F8', zorder=0)
    ax.add_feature(cfeature.LAND,  facecolor='#FAFAF5', zorder=0)
    if title:
        ax.set_title(title, fontsize=9, fontweight='bold', pad=4)

def add_features(ax):
    ax.add_feature(cfeature.OCEAN, facecolor='#E8F4F8', zorder=10)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor='black', zorder=11)
    ax.add_feature(MA_BORDERS, linewidth=0.4, edgecolor='black', zorder=11)

# ============================================================================
def extract_score_data(ds, score_name, variable, category_idx=None):
    """
    Extraire les donnees d'un score depuis un dataset OSOP.
    - variable : 't2m' ou 'tp'
    - category_idx : 0=below, 1=normal, 2=above (None si pas de category)
    Retour : 2D array (lat, lon)
    """
    var_name = 't2m' if variable == '2m_temperature' else 'tp'
    if var_name not in ds.data_vars:
        var_name = list(ds.data_vars)[0]  # fallback

    arr = ds[var_name]

    # Squeeze forecastMonth
    if 'forecastMonth' in arr.dims:
        arr = arr.isel(forecastMonth=0)

    # Selectionner category si applicable
    if 'category' in arr.dims:
        if category_idx is None:
            category_idx = 2  # default = above
        arr = arr.isel(category=category_idx)

    return arr.values

# ============================================================================
# CONFIGURATION DES SCORES A GENERER
# ============================================================================
score_configs = [
    {
        'score': 'spearman_corr',
        'category_idx': None,
        'label_short': 'Spearman',
        'label_full': 'Spearman Rank Correlation',
        'cmap': plt.cm.RdBu_r,
        'vmin': -1.0, 'vmax': 1.0,
        'cbar_extend': 'neither',
    },
    {
        'score': 'bs',
        'category_idx': 2,  # above-normal
        'label_short': 'BS_above',
        'label_full': 'Brier Score (above-normal) — lower is better',
        'cmap': plt.cm.RdYlGn_r,
        'vmin': 0.0, 'vmax': 0.4,
        'cbar_extend': 'max',
    },
    {
        'score': 'rps',
        'category_idx': None,
        'label_short': 'RPS',
        'label_full': 'Ranked Probability Score — lower is better',
        'cmap': plt.cm.RdYlGn_r,
        'vmin': 0.0, 'vmax': 1.5,
        'cbar_extend': 'max',
    },
    {
        'score': 'roc',
        'category_idx': 2,  # above-normal
        'label_short': 'ROC_above',
        'label_full': 'ROC Area (above-normal) — 0.5 = no skill',
        'cmap': plt.cm.RdYlGn,
        'vmin': 0.3, 'vmax': 1.0,
        'cbar_extend': 'min',
    },
    {
        'score': 'rocss',
        'category_idx': 0,  # below-normal
        'label_short': 'ROCSS_below',
        'label_full': 'ROC Skill Score (below-normal)',
        'cmap': plt.cm.RdYlGn,
        'vmin': -0.5, 'vmax': 0.5,
        'cbar_extend': 'both',
    },
]

# ============================================================================
# CONFIGURATION DES VARIABLES
# ============================================================================
variable_configs = [
    {
        'variable':  'total_precipitation',
        'var_short': 'TP',
        'subdir':    'Precipitation',
        'fig_prefix':'Fig5bis_OSOP_multisystem_skill_TP',
    },
    {
        'variable':  '2m_temperature',
        'var_short': 'T2m',
        'subdir':    'TMP',
        'fig_prefix':'Fig5ter_OSOP_multisystem_skill_T2m',
    },
]

# ============================================================================
# BOUCLE PRINCIPALE
# ============================================================================
total_figures = len(score_configs) * len(variable_configs)
fig_count = 0

for var_cfg in variable_configs:
    for sc_cfg in score_configs:
        fig_count += 1
        print(f"\n{'='*70}")
        print(f"[{fig_count}/{total_figures}] {var_cfg['var_short']} - {sc_cfg['label_short']}")
        print(f"{'='*70}")

        # Charger les datasets
        datasets = {}
        for prefix in MODEL_ORDER:
            f = find_score_file(prefix, var_cfg['variable'], sc_cfg['score'])
            if f is None:
                print(f"  [MISSING] {prefix:18s}")
                datasets[prefix] = None
            else:
                try:
                    datasets[prefix] = xr.open_dataset(f)
                    print(f"  [OK]      {prefix:18s} -> {os.path.basename(f)}")
                except Exception as e:
                    print(f"  [ERROR]   {prefix:18s} -> {str(e)[:60]}")
                    datasets[prefix] = None

        n_found = sum(1 for d in datasets.values() if d is not None)
        if n_found == 0:
            print(f"  Skip : aucun fichier trouve.")
            continue

        # Generer la figure
        fig = plt.figure(figsize=(18, 9))
        gs = GridSpec(2, 5, figure=fig, hspace=0.28, wspace=0.10,
                      left=0.03, right=0.94, top=0.91, bottom=0.10)

        last_im = None
        for idx, prefix in enumerate(MODEL_ORDER):
            row, col = idx // 5, idx % 5
            ax = fig.add_subplot(gs[row, col], projection=ccrs.PlateCarree())
            label = MODEL_LABELS.get(prefix, prefix.upper())
            setup_subplot(ax, title=label)

            ds = datasets.get(prefix)
            if ds is not None:
                try:
                    data_2d = extract_score_data(ds, sc_cfg['score'], var_cfg['variable'],
                                                  category_idx=sc_cfg['category_idx'])
                    # Lat/lon
                    lat_name = 'lat'
                    lon_name = 'lon'
                    im = ax.pcolormesh(ds[lon_name], ds[lat_name], data_2d,
                                        transform=ccrs.PlateCarree(),
                                        cmap=sc_cfg['cmap'],
                                        vmin=sc_cfg['vmin'], vmax=sc_cfg['vmax'],
                                        shading='auto', zorder=2)
                    last_im = im
                    add_features(ax)
                except Exception as e:
                    ax.text(0.5, 0.5, f"Error\n{str(e)[:30]}",
                            ha='center', va='center', transform=ax.transAxes,
                            fontsize=8, color='red', style='italic')
            else:
                ax.text(0.5, 0.5, f"File not\nfound",
                        ha='center', va='center', transform=ax.transAxes,
                        fontsize=10, color='gray', style='italic')

        # Colorbar verticale a droite
        if last_im is not None:
            cax = fig.add_axes([0.955, 0.15, 0.013, 0.65])
            cb = plt.colorbar(last_im, cax=cax, orientation='vertical',
                              extend=sc_cfg['cbar_extend'])
            cb.set_label(sc_cfg['label_full'], fontsize=9, fontweight='bold')
            cb.ax.tick_params(labelsize=8)

        fig.suptitle(f"Multi-System Skill — {var_cfg['var_short']} JJA Hindcast 1993-2016 — North Africa\n"
                      f"{sc_cfg['label_full']}",
                      fontsize=14, fontweight='bold', y=0.97)

        fig.text(0.5, 0.025,
                  f"C3S seasonal forecast systems (OSOP/SEAFOAM pipeline) — Hindcast period: 1993-2016  |  "
                  f"Verification: ERA5  |  "
                  f"Note: JMA CPS4 shown for reference but excluded from the MME.",
                  ha='center', fontsize=8.5, style='italic', color='#555555')

        # Sauvegarder
        fname = f"{var_cfg['fig_prefix']}_{sc_cfg['label_short']}_JJA.png"
        out_file = os.path.join(OUT_DIR, fname)
        plt.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.2)

        # Copie dans JJA2026/<subdir>/
        out_file_jja = os.path.join(JJA_DIR, var_cfg['subdir'], fname)
        shutil.copy(out_file, out_file_jja)

        plt.close()
        print(f"\n[OK] {fname}")
        print(f"     {os.path.getsize(out_file)/1024:.0f} KB")
        print(f"[OK] Copie aussi dans : {out_file_jja}")

        # Fermer datasets
        for ds in datasets.values():
            if ds is not None:
                ds.close()

print("\n" + "=" * 70)
print(f"  GENERATION TERMINEE : {total_figures} figures multi-modele skill")
print(f"  Voir : {JJA_DIR}")
