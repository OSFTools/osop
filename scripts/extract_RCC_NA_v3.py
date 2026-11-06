#!/usr/bin/env python
"""
Extraction des chiffres clés pour le bulletin RCC Afrique du Nord — JJA 2026.

Cette version produit DEUX vues :
  1. Par CAPITALE (point géographique)
  2. Par PAYS (moyenne sur bounding box)

Sortie : CSV consolidés prêts pour insertion dans le Word.
"""
import glob
import os
import numpy as np
import pandas as pd
import xarray as xr

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE = 'C:/Users/wafae/projets/osop-data/seafoam/data/master'
FORECAST_DIR = f'{BASE}/forecast/products'
HINDCAST_SCORES_DIR = f'{BASE}/hindcast/scores'
OUTPUT_CSV_FORECAST = f'{BASE}/RCC_NA_JJA2026_forecast.csv'
OUTPUT_CSV_SKILL    = f'{BASE}/RCC_NA_JJA2026_skill.csv'
OUTPUT_CSV_COUNTRY  = f'{BASE}/RCC_NA_JJA2026_by_country.csv'

# ----------------------------------------------------------------------------
# 1. Capitales (points)
# ----------------------------------------------------------------------------
CITIES = {
    'Rabat (MA)':       (34.0,  -6.8),
    'Casablanca (MA)':  (33.6,  -7.6),
    'Marrakech (MA)':   (31.6,  -8.0),
    'Oujda (MA)':       (34.7,  -1.9),
    'Alger (DZ)':       (36.7,   3.1),
    'Oran (DZ)':        (35.7,  -0.6),
    'Constantine (DZ)': (36.4,   6.6),
    'Tamanrasset (DZ)': (22.8,   5.5),
    'Tunis (TN)':       (36.8,  10.2),
    'Sfax (TN)':        (34.7,  10.8),
    'Tripoli (LY)':     (32.9,  13.2),
    'Benghazi (LY)':    (32.1,  20.1),
    'Le Caire (EG)':    (30.0,  31.2),
    'Alexandrie (EG)':  (31.2,  29.9),
    'Nouakchott (MR)':  (18.1, -16.0),
    'Khartoum (SD)':    (15.6,  32.5),
}

# ----------------------------------------------------------------------------
# 2. Bounding boxes pays (lat_min, lat_max, lon_min, lon_max)
# Sources : limites terrestres approximatives — restent dans le domaine OSOP.
# ----------------------------------------------------------------------------
COUNTRY_BOXES = {
    'Morocco':       {'lat': (27.5, 35.9), 'lon': (-13.0, -1.0)},
    'Algeria North': {'lat': (32.0, 37.0), 'lon': ( -8.0, 12.0)},  # Tell + Hauts Plateaux
    'Algeria South': {'lat': (19.0, 28.0), 'lon': ( -5.0, 11.0)},  # Sahara
    'Tunisia':       {'lat': (30.0, 37.4), 'lon': (  7.5, 11.6)},
    'Libya':         {'lat': (19.5, 33.2), 'lon': (  9.4, 25.2)},
    'Egypt':         {'lat': (22.0, 31.7), 'lon': ( 24.7, 36.9)},
    'Mauritania':    {'lat': (14.7, 27.3), 'lon': (-17.1, -4.8)},
    'Maghreb':       {'lat': (28.0, 37.0), 'lon': (-13.0, 12.0)},  # synthèse régionale
    'NorthAfrica':   {'lat': (15.0, 37.0), 'lon': (-17.0, 37.0)},  # tout le RCC
}

# Centres : MME en priorité, puis quelques majors pour convergence
CENTRES_FORECAST = {
    'mme': 'MME multi-modèle',
    'ecmwf_51': 'ECMWF SEAS5',
    'ukmo_604': 'UKMO GloSea6',
    'meteo_france_9': 'Météo-France 9',
}

# ============================================================================
# UTILS
# ============================================================================
def find_file(pattern, directory):
    files = glob.glob(os.path.join(directory, pattern))
    return files[0] if files else None

def sel_point(ds, lat, lon):
    lat_name = 'latitude' if 'latitude' in ds.dims else 'lat'
    lon_name = 'longitude' if 'longitude' in ds.dims else 'lon'
    return ds.sel({lat_name: lat, lon_name: lon}, method='nearest')

def sel_box(ds, lat_min, lat_max, lon_min, lon_max):
    """Sélectionne une boîte lat/lon, gère le sens du lat (croissant ou décroissant)."""
    lat_name = 'latitude' if 'latitude' in ds.dims else 'lat'
    lon_name = 'longitude' if 'longitude' in ds.dims else 'lon'
    lat_vals = ds[lat_name].values
    if lat_vals[0] > lat_vals[-1]:  # décroissant
        slice_lat = slice(lat_max, lat_min)
    else:
        slice_lat = slice(lat_min, lat_max)
    return ds.sel({lat_name: slice_lat, lon_name: slice(lon_min, lon_max)})

def mean_box(ds, lat_min, lat_max, lon_min, lon_max, var):
    sub = sel_box(ds, lat_min, lat_max, lon_min, lon_max)
    return float(sub[var].mean().values)

# ============================================================================
# 1. FORECAST PAR CAPITALE
# ============================================================================
print("\n" + "=" * 90)
print("1. PROBABILITÉS TERCILES — PAR CAPITALE")
print("=" * 90)

forecast_results = []
for centre_id, centre_label in CENTRES_FORECAST.items():
    for var_name, var_label in [('2m_temperature', 'T2m'), ('total_precipitation', 'TP')]:
        f = find_file(f'{centre_id}_*{var_name}.3m.forecast_percentages.nc', FORECAST_DIR)
        if not f:
            continue
        ds = xr.open_dataset(f)
        for city, (lat, lon) in CITIES.items():
            try:
                pt = sel_point(ds, lat, lon)
                p_low  = float(pt['lower'].values.flatten()[0])
                p_mid  = float(pt['middle'].values.flatten()[0])
                p_high = float(pt['higher'].values.flatten()[0])
                clim = 33.3
                if p_high - clim > 7 and p_high > p_low:
                    signal = 'ABOVE'
                elif p_low - clim > 7 and p_low > p_high:
                    signal = 'BELOW'
                else:
                    signal = 'normal/incertain'
                forecast_results.append({
                    'Centre': centre_label, 'Variable': var_label, 'Ville': city,
                    'P_below(%)':  round(p_low, 1),
                    'P_normal(%)': round(p_mid, 1),
                    'P_above(%)':  round(p_high,1),
                    'Signal': signal,
                })
            except Exception as e:
                print(f"   ⚠️  {centre_id} {var_label} {city} : {e}")
        ds.close()
df_forecast = pd.DataFrame(forecast_results)

df_mme_cap = df_forecast[df_forecast['Centre'] == 'MME multi-modèle'].copy()
print("\n📊 MME par capitale :")
print(df_mme_cap[['Variable','Ville','P_below(%)','P_normal(%)','P_above(%)','Signal']]
      .to_string(index=False))

# ============================================================================
# 2. FORECAST PAR PAYS (MOYENNES DE BOX)
# ============================================================================
print("\n" + "=" * 90)
print("2. PROBABILITÉS TERCILES — MOYENNE PAR PAYS / RÉGION")
print("=" * 90)

country_results = []
for centre_id, centre_label in CENTRES_FORECAST.items():
    for var_name, var_label in [('2m_temperature', 'T2m'), ('total_precipitation', 'TP')]:
        f = find_file(f'{centre_id}_*{var_name}.3m.forecast_percentages.nc', FORECAST_DIR)
        if not f:
            continue
        ds = xr.open_dataset(f)
        for country, box in COUNTRY_BOXES.items():
            try:
                lat_min, lat_max = box['lat']
                lon_min, lon_max = box['lon']
                p_low  = mean_box(ds, lat_min, lat_max, lon_min, lon_max, 'lower')
                p_mid  = mean_box(ds, lat_min, lat_max, lon_min, lon_max, 'middle')
                p_high = mean_box(ds, lat_min, lat_max, lon_min, lon_max, 'higher')
                clim = 33.3
                if p_high - clim > 5 and p_high > p_low:
                    signal = 'ABOVE'
                elif p_low - clim > 5 and p_low > p_high:
                    signal = 'BELOW'
                else:
                    signal = 'normal/incertain'
                country_results.append({
                    'Centre': centre_label, 'Variable': var_label, 'Pays': country,
                    'P_below(%)':  round(p_low, 1),
                    'P_normal(%)': round(p_mid, 1),
                    'P_above(%)':  round(p_high,1),
                    'Signal': signal,
                })
            except Exception as e:
                print(f"   ⚠️  {centre_id} {var_label} {country} : {e}")
        ds.close()
df_country = pd.DataFrame(country_results)

df_mme_country = df_country[df_country['Centre'] == 'MME multi-modèle'].copy()
print("\n📊 MME moyenne par pays :")
print(df_mme_country[['Variable','Pays','P_below(%)','P_normal(%)','P_above(%)','Signal']]
      .to_string(index=False))

# ============================================================================
# 3. SKILL HINDCAST — MOYENNE PAR PAYS
# ============================================================================
print("\n" + "=" * 90)
print("3. SKILL HINDCAST — MOYENNE PAR PAYS")
print("=" * 90)

skill_country = []
for var_name, var_label in [('2m_temperature', 'T2m'), ('total_precipitation', 'TP')]:
    fp = find_file(f'mme*{var_name}*3m*pearson_corr.nc', HINDCAST_SCORES_DIR)
    fr = find_file(f'mme*{var_name}*3m*rocss.nc',        HINDCAST_SCORES_DIR)
    if not fp:
        continue
    dsp = xr.open_dataset(fp)
    dsr = xr.open_dataset(fr) if fr else None
    var_p = list(dsp.data_vars)[0]
    var_r = list(dsr.data_vars)[0] if dsr else None

    for country, box in COUNTRY_BOXES.items():
        try:
            lat_min, lat_max = box['lat']
            lon_min, lon_max = box['lon']
            pearson = mean_box(dsp, lat_min, lat_max, lon_min, lon_max, var_p)
            rocss_b, rocss_a = None, None
            if dsr is not None:
                sub = sel_box(dsr, lat_min, lat_max, lon_min, lon_max)
                arr = sub[var_r].values
                # arr shape attendue : (category, lat, lon)  ou  (lat, lon, category)
                if arr.ndim == 3 and arr.shape[0] == 3:
                    rocss_b = float(np.nanmean(arr[0]))
                    rocss_a = float(np.nanmean(arr[2]))
                elif arr.ndim == 3 and arr.shape[-1] == 3:
                    rocss_b = float(np.nanmean(arr[..., 0]))
                    rocss_a = float(np.nanmean(arr[..., 2]))
            niveau = ('🟢 Robuste' if pearson > 0.5
                      else '🟡 Modérée' if pearson > 0.3
                      else '🔴 Faible')
            skill_country.append({
                'Variable': var_label, 'Pays': country,
                'Pearson_r':    round(pearson, 2),
                'ROCSS_below':  round(rocss_b, 2) if rocss_b is not None else None,
                'ROCSS_above':  round(rocss_a, 2) if rocss_a is not None else None,
                'Niveau_skill': niveau,
            })
        except Exception as e:
            print(f"   ⚠️  {country} {var_label} : {e}")
    dsp.close()
    if dsr is not None: dsr.close()

df_skill_country = pd.DataFrame(skill_country)
print("\n📈 Skill MME moyenne par pays :")
print(df_skill_country.to_string(index=False))

# ============================================================================
# 4. EXPORTS
# ============================================================================
print("\n" + "=" * 90)
print("EXPORT")
print("=" * 90)
df_forecast.to_csv(OUTPUT_CSV_FORECAST, index=False, encoding='utf-8-sig')
df_country.to_csv(OUTPUT_CSV_COUNTRY, index=False, encoding='utf-8-sig')
df_skill_country.to_csv(OUTPUT_CSV_SKILL.replace('.csv','_country.csv'),
                        index=False, encoding='utf-8-sig')
print(f"✅ {OUTPUT_CSV_FORECAST}        ({len(df_forecast)} lignes)")
print(f"✅ {OUTPUT_CSV_COUNTRY}          ({len(df_country)} lignes)")
print(f"✅ {OUTPUT_CSV_SKILL.replace('.csv','_country.csv')}   ({len(df_skill_country)} lignes)")
print("\n✅ EXTRACTION TERMINÉE\n")
