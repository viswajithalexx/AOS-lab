#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 19:19:42 2026

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


file = '/home/bobco-08/24cl05012/CO2/data/data_1/cmems/bio_var/cmems_obs-mob_glo_bgc-car_my_irr-i_multi-vars_35.12E-119.88E_29.88S-29.88N_1994-01-01-2024-12-01.nc'

data = xr.open_dataset(file)
data = data.rename({'latitude':'lat', 'longitude': 'lon'})
lat = data['lat']
lon =data['lon']
time = data['time']
data2 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/cmems/phy_var/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_multi-vars_35.00E-120.00E_30.00S-30.00N_0.51m_1993-01-01-2024-12-01.nc')
data2 = data2.rename({'latitude':'lat', 'longitude': 'lon'})
sst = data2['thetao_oras'].isel(depth=0)
sss = data2['so_oras'].isel(depth=0)



#%%
regions = {
    'NWIO': dict(lat=slice(5,22.5), lon=slice(45,65)),
    'ESIO': dict(lat=slice(-6.6,8), lon=slice(92,109)),
    'EIO':  dict(lat=slice(-6.5,5), lon=slice(49,92)),
    'NAS':  dict(lat=slice(22.5,28), lon=slice(56,70)),
}
def compute_pco2_terms(region_name, lat_slice, lon_slice):

    # ---- Extract region ----
    dic = data['tco2'].sel(lat=lat_slice, lon=lon_slice)
    alk = data['talk'].sel(lat=lat_slice, lon=lon_slice)
    pco2 = data['spco2'].sel(lat=lat_slice, lon=lon_slice)

    sst_r = sst.sel(lat=lat_slice, lon=lon_slice,
                    time=slice('1994-01-01','2024-12-01'))
    sss_r = sss.sel(lat=lat_slice, lon=lon_slice,
                    time=slice('1994-01-01','2024-12-01'))

    # ---- Mean over region ----
    pco2_reg = pco2.mean(dim=('lat','lon'))
    dic_reg  = dic.mean(dim=('lat','lon'))
    alk_reg  = alk.mean(dim=('lat','lon'))
    sst_reg  = sst_r.mean(dim=('lat','lon'))
    sal_reg  = sss_r.mean(dim=('lat','lon'))

    # ---- Annual ----
    pco2_annual = pco2_reg.groupby('time.year').mean()
    dic_annual  = dic_reg.groupby('time.year').mean()
    alk_annual  = alk_reg.groupby('time.year').mean()

    # ---- Seasonal helper ----
    def seasonal_mean(var, months):
        return var.sel(time=var['time'].dt.month.isin(months)).groupby('time.year').mean()

    # ---- JJA ----
    pco2_jja = seasonal_mean(pco2_reg, [7,8,9])
    dic_jja  = seasonal_mean(dic_reg, [7,8,9])
    alk_jja  = seasonal_mean(alk_reg, [7,8,9])
    sal_jja  = seasonal_mean(sal_reg, [7,8,9])
    sst_jja  = seasonal_mean(sst_reg, [7,8,9])

    # ---- DJF ----
    def djf_mean(var):
        v = var.resample(time='QS-DEC').mean()
        v = v.isel(time=slice(1,-1))
        v = v.sel(time=v.time.dt.month == 12)
        v['year'] = v['time.year'] + 1
        return v.groupby('year').mean()

    pco2_djf = djf_mean(pco2_reg)
    dic_djf  = djf_mean(dic_reg)
    alk_djf  = djf_mean(alk_reg)
    sal_djf  = djf_mean(sal_reg)
    sst_djf  = djf_mean(sst_reg)

    # ---- Differences ----
    delta_pco2 = pco2_jja - pco2_djf
    delta_T    = sst_jja - sst_djf
    delta_dic  = dic_jja - dic_djf
    delta_alk  = alk_jja - alk_djf
    delta_sal  = sal_jja - sal_djf

    # ---- Terms ----
    T_term   = 2 * pco2_annual * (np.exp(0.0423 * (delta_T / 2)) - 1)
    dic_term = (9.5 * (pco2_annual/dic_annual)) * delta_dic
    alk_term = (-8.9 * (pco2_annual/alk_annual)) * delta_alk
    sal_term = (0.026 * pco2_annual) * delta_sal

    delta_cal = T_term + dic_term + alk_term + sal_term

    # ---- Align ----
    delta_pco2, delta_cal = xr.align(delta_pco2, delta_cal)

    # ---- RMSE ----
    rmse = np.sqrt(np.mean((delta_cal - delta_pco2)**2))
    print(f"{region_name} RMSE = {rmse.values:.2f} µatm")

    # =========================
    # MEAN CONTRIBUTIONS
    # =========================
    mean_T   = float(T_term.mean().values)
    mean_DIC = float(dic_term.mean().values)
    mean_ALK = float(alk_term.mean().values)
    mean_SAL = float(sal_term.mean().values)
    mean_CAL = float(delta_cal.mean().values)
    mean_OBS = float(delta_pco2.mean().values)

    total = mean_CAL if mean_CAL != 0 else np.nan

    contrib_T   = (mean_T / total) * 100
    contrib_DIC = (mean_DIC / total) * 100
    contrib_ALK = (mean_ALK / total) * 100
    contrib_SAL = (mean_SAL / total) * 100

    # =========================
    # PRINT SUMMARY
    # =========================
    print(f"\n===== {region_name} =====")

    print("Mean Contributions (µatm):")
    print(f"T   : {mean_T:.2f}")
    print(f"DIC : {mean_DIC:.2f}")
    print(f"ALK : {mean_ALK:.2f}")
    print(f"SAL : {mean_SAL:.2f}")
    print(f"CAL : {mean_CAL:.2f}")
    print(f"OBS : {mean_OBS:.2f}")

    print("\nContribution (%):")
    print(f"T   : {contrib_T:.1f}%")
    print(f"DIC : {contrib_DIC:.1f}%")
    print(f"ALK : {contrib_ALK:.1f}%")
    print(f"SAL : {contrib_SAL:.1f}%")

    print("\nDifference (CAL - OBS):")
    print(f"{(mean_CAL - mean_OBS):.2f} µatm")

    # ---- RETURN ----
    return {
        'years': delta_pco2['year'].values,
        'obs': delta_pco2.values,
        'cal': delta_cal.values,
        'T': T_term.values,
        'DIC': dic_term.values,
        'ALK': alk_term.values,
        'SAL': sal_term.values
    }
results = {}

for name, reg in regions.items():
    results[name] = compute_pco2_terms(name, reg['lat'], reg['lon'])
    
#%%
import numpy as np
import matplotlib.pyplot as plt

# =========================
# PREPARE SUMMARY (USE EXISTING results)
# =========================
summary = {}

for name, res in results.items():

    T   = res['T']
    DIC = res['DIC']
    ALK = res['ALK']
    SAL = res['SAL']
    CAL = res['cal']
    OBS = res['obs']

    means = np.array([
        np.mean(T),
        np.mean(DIC),
        np.mean(ALK),
        np.mean(SAL),
        np.mean(CAL),
        np.mean(OBS)
    ])

    stds = np.array([
        np.std(T, ddof=1),
        np.std(DIC, ddof=1),
        np.std(ALK, ddof=1),
        np.std(SAL, ddof=1),
        np.std(CAL, ddof=1),
        np.std(OBS, ddof=1)
    ])

    total = np.mean(CAL)
    contrib = [(m/total)*100 if total != 0 else 0 for m in means[:4]]

    summary[name] = {
        'means': means,
        'stds': stds,
        'contrib': contrib
    }

# =========================
# PLOTTING
# =========================
variables = ['T', 'DIC', 'ALK', 'SAL', 'CAL', 'OBS']
colors = ['red', 'green', 'purple', 'blue', 'grey', 'lightblue']

fig, axes = plt.subplots(2,2, figsize=(20,18), dpi=600)
axes = axes.flatten()

for idx, (region, data) in enumerate(summary.items()):

    ax = axes[idx]

    means = data['means']
    errors = data['stds']
    contrib = data['contrib']

    x = np.arange(len(variables))

    bars = []

    # ---- BAR PLOT ----
    for i in range(len(variables)):
        bar = ax.bar(
        x[i], means[i],
        color=colors[i],
        width=0.5,
        yerr=errors[i],
        capsize=8,
        ecolor='black',
        edgecolor='black',   # ✅ adds black boundary
        linewidth=1.5        # ✅ thickness of border
    )

    # ---- ZERO LINE ----
    ax.axhline(0, color='black', linewidth=1)

    # ---- TITLE ----
    ax.set_title(region, fontsize=24, fontweight='bold')

    # ---- Y LABEL ----
    # if idx % 2 == 0:
    ax.set_ylabel('ΔpCO$_2$ (µatm)', fontsize=24)
    ax.tick_params(axis='y', labelsize=24)
    # ---- X LABELS ----
    ax.set_xticks(x)
    ax.set_xticklabels(variables, fontsize=24, fontweight='bold')

    # ---- GRID ----
    ax.grid(axis='y', linestyle='--', alpha=0.3)

    # =========================
    # VALUES near x-axis
    # =========================
    y_offset = -0.25 * np.max(np.abs(means))

    # for i, bar in enumerate(bars):
    #     val = means[i]

        # ax.text(
        #     bar.get_x() + bar.get_width()/2,
        #     y_offset,
        #     f'{val:.2f}',
        #     ha='center',
        #     fontsize=11,
        #     fontweight='bold',
        #     color='red' if val < 0 else 'black'
        # )
# ---- REMOVE EXTRA AXES ----
for j in range(idx+1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.show()