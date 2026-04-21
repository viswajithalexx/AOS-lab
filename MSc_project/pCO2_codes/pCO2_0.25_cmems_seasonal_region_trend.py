#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 18:23:50 2026

@author: bobco-08
"""

import xarray as xr
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# =========================
# 1. Load Data
# =========================
file1 = "/home/bobco-08/24cl05012/CO2/data/data_1/jena/oc_v2025.flux.nc"

ds = xr.open_dataset(file1)
# ds = ds.rename({'latitude': 'lat', 'longitude': 'lon'})

ds = ds.rename({'mtime': 'time'})
time = ds['time']
dxyp = ds['dxyp']
co2f = ds['co2flux_ocean']
co2f = co2f/dxyp

co2fv = co2f.values*1e15
# for IO

pco2 = co2f*1e15
#%%
# =========================
# 2. Define Regions
# =========================
regions = { 
    'NWIO':  dict(lat=slice(5,22.5), lon=slice(45,65)),         #selecting regions 
    'NAS':  dict(lat=slice(22.5,28), lon=slice(56,70)),
    'EIO': dict(lat = slice(-6.5,5),lon = slice(49,92)),
    'ESIO': dict(lat=slice(-6.6,8), lon=slice(92,109)),
}
#%%
# =========================
# 3. Seasonal Function
# =========================
def seasonal_timeseries(data, months, season_name):
    
    if season_name == 'DJF':                                   #selecting DJF
        ds_season = data.resample(time='QS-DEC').mean()
        ds_season = ds_season.sel(time=ds_season.time.dt.month == 12)
    else:
        ds_season = data.sel(time=data.time.dt.month.isin(months)) #selecting JJAS,MAM,ON
        ds_season = ds_season.resample(time='YE').mean()

    # spatial mean
    ds_mean = ds_season.mean(dim=('lat', 'lon'))

    # remove NaNs
    ds_mean = ds_mean.dropna(dim='time')

    return ds_mean
#%%
# =========================
# 4. Seasons
# =========================
seasons = {
    'DJF':  [12,1,2],
    'MAM':  [3,4,5],
    'JJAS': [6,7,8,9],
    'ON':   [10,11]
}
#%%
def plot_region(region_name, region_data, ylim=None):

    fig, axs = plt.subplots(2, 2, figsize=(16, 12))
    axs = axs.flatten()

    for i, (season, months) in enumerate(seasons.items()):

        ax = axs[i]

        ts = seasonal_timeseries(region_data, months, season)

        years = ts.time.dt.year.values
        values = ts.values

        # Trend
        coeff = np.polyfit(years, values, 1)
        trend = np.polyval(coeff, years)

        stats_out = stats.linregress(years, values)
        pval = stats_out.pvalue

        label = f"({coeff[0]:.2f} µatm/yr,p={pval:.3f})"

        # Plot
        ax.plot(years, values, 'o-')
        ax.plot(years, trend, 'r--', lw=2, label=label)

        # ===== TITLES =====
        ax.set_title(season, fontsize=22, fontweight='bold')

        # ===== AXIS LABELS =====
        if i in [0, 2]:
            ax.set_ylabel('µatm', fontsize=20, fontweight='bold')

        if i in [2, 3]:
            ax.set_xlabel('Years', fontsize=20, fontweight='bold')

        # ===== Y LIMIT =====
        if ylim is not None:
            ax.set_ylim(ylim)

        # ===== STYLE =====
        ax.tick_params(axis='both', labelsize=16)
        ax.grid(True, linestyle=':')

        # Legend in all panels
        ax.legend(fontsize=14)

    # ===== FIGURE TITLE =====
    fig.suptitle(f"{region_name} pCO$_2$ Trends (1957–2024)", 
                 fontsize=24, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()
    
    
region_data_dict = {}

for name, bounds in regions.items():

    region_data = pco2.sel(**bounds)

    plot_region(name, region_data)
    
#%% 

import xarray as xr
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# =========================
# 1. Load Data
# =========================
file1 = "/home/bobco-08/24cl05012/CO2/data/data_1/jena/oc_v2025.pCO2.nc"

ds = xr.open_dataset(file1)
# ds = ds.rename({'latitude': 'lat', 'longitude': 'lon'})
ds = ds.rename({'mtime': 'time'})
pco2= ds['pCO2'] 

# Time selection
pco2 = pco2.sel(time=slice('1957-01-01', '2024-12-31'))

# =========================
# 2. Define Regions
# =========================
regions = { 
    'NWIO':  dict(lat=slice(5,22.5), lon=slice(45,65)),         #selecting regions 
    'NAS':  dict(lat=slice(22.5,28), lon=slice(56,70)),
    'EIO': dict(lat = slice(-6.5,5),lon = slice(49,92)),
    'ESIO': dict(lat=slice(-6.6,8), lon=slice(92,109)),
}
#%%
# =========================
# 3. Seasonal Function
# =========================
def seasonal_timeseries(data, months, season_name):
    
    if season_name == 'DJF':                                   #selecting DJF
        ds_season = data.resample(time='QS-DEC').mean()
        ds_season = ds_season.sel(time=ds_season.time.dt.month == 12)
    else:
        ds_season = data.sel(time=data.time.dt.month.isin(months)) #selecting JJAS,MAM,ON
        ds_season = ds_season.resample(time='YE').mean()

    # spatial mean
    ds_mean = ds_season.mean(dim=('lat', 'lon'))

    # remove NaNs
    ds_mean = ds_mean.dropna(dim='time')

    return ds_mean
#%%
# =========================
# 4. Seasons
# =========================
seasons = {
    'DJF':  [12,1,2],
    'MAM':  [3,4,5],
    'JJAS': [6,7,8,9],
    'ON':   [10,11]
}
#%%
# =========================
# 5. Plot Function
# =========================
def plot_all_regions(regions_dict, ylim=None):

    fig, axs = plt.subplots(4, 4, figsize=(21,18))

    for r, (region_name, region_data) in enumerate(regions_dict.items()):

        for c, (season, months) in enumerate(seasons.items()):

            ax = axs[r, c]

            ts = seasonal_timeseries(region_data, months, season)

            years = ts.time.dt.year.values
            values = ts.values

            # Trend
            coeff = np.polyfit(years, values, 1)
            trend = np.polyval(coeff, years)

            stats_out = stats.linregress(years, values)
            pval = stats_out.pvalue

            label = f"({coeff[0]:.2f} µatm/yr)"

            # Plot
            ax.plot(years, values, 'o-')
            ax.plot(years, trend, 'r--', lw=2, label=label)

            # ===== TITLES =====
            if r == 0:
                ax.set_title(season, fontsize=21, fontweight='bold')

            if c == 0:
                ax.set_ylabel(f"{region_name}\nµatm", fontsize=21, fontweight = 'bold')

            # ===== AXIS LABELS =====
            if r == 3:
                ax.set_xlabel('Years', fontsize=21,fontweight = 'bold')
            # ===== Y LIMIT OPTION =====
            if ylim is not None:
                ax.set_ylim(ylim)

            # ===== STYLE =====
            ax.tick_params(axis='both', labelsize=21)
            ax.grid(True, linestyle=':')

            # Legend only in last column
            ax.legend(fontsize=20)
    # ===== LAYOUT =====
    # fig.suptitle("Seasonal pCO$_2$ Trends (1957–2024)", fontsize=28, y=0.98)

    plt.subplots_adjust(
        left=0.08,
        right=0.97,
        top=0.93,
        bottom=0.08,
        wspace=0.19,
        hspace=0.22
    )

    plt.show()
    
    
region_data_dict = {}

for name, bounds in regions.items():
    region_data_dict[name] = pco2.sel(**bounds)

plot_all_regions(region_data_dict, ylim=(220,450))
#%%
def plot_region(region_name, region_data, ylims=None):

    fig, axs = plt.subplots(2, 2, figsize=(20, 10))
    axs = axs.flatten()

    for i, (season, months) in enumerate(seasons.items()):

        ax = axs[i]

        ts = seasonal_timeseries(region_data, months, season)

        years = ts.time.dt.year.values
        values = ts.values

        # Trend
        coeff = np.polyfit(years, values, 1)
        trend = np.polyval(coeff, years)

        stats_out = stats.linregress(years, values)
        pval = stats_out.pvalue

        label = f"({coeff[0]*10:.2f} gC m⁻² yr⁻¹ decade⁻¹), p={pval:.3f})"

        # Plot
        ax.plot(years, values, 'o-', lw=3.5)
        ax.plot(years, trend, 'r--', lw=2, label=label)

        # ===== TITLES =====
        ax.set_title(f"{region_name} - {season}", fontsize=24, pad=15, fontweight ='bold')

        # ===== AXIS LABELS =====
        if i in [0, 2]:
            ax.set_ylabel('flux (gC m$^{-2}$ yr$^{-1}$)', fontsize=24,labelpad=15)

        if i in [2, 3]:
            ax.set_xlabel('Years', fontsize=24,labelpad=15)

        # ===== MANUAL Y-LIMIT PER SEASON =====
        if ylims is not None and season in ylims:
            ax.set_ylim(ylims[season])

        # ===== STYLE =====
        ax.tick_params(axis='both', labelsize=24)
        ax.grid(True, linestyle=':')
        for side in ['top', 'bottom', 'left', 'right']:
            ax.spines[side].set_linewidth(2.5)
        ax.legend(fontsize=24)

    # ===== LAYOUT =====
    plt.subplots_adjust(
        left=0.08,
        right=0.97,
        top=0.97,
        bottom=0.08,
        wspace=0.11,
        hspace=0.27
    )
    plt.show()



ylims_all = {
    'NWIO': {'DJF': (-4,12), 'MAM': (-4,12), 'JJAS':None, 'ON': (-4,12)},
    'NAS':  {'DJF': (0,25), 'MAM': (0,25), 'JJAS': (0,25), 'ON': (0,25)},
    'EIO':  {'DJF': (0,12), 'MAM': (0,12), 'JJAS': (0,12), 'ON': (0,12)},
    'ESIO': {'DJF': (-6,6), 'MAM': (-6,6), 'JJAS': (-6,6), 'ON': (-6,6)},
}

for name, bounds in regions.items():
    region_data = pco2.sel(**bounds)

    plot_region(name, region_data, ylims=ylims_all[name])
