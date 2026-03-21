#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 12:34:07 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


file1 = "/home/bobco-08/24cl05012/CO2/data/oc_v2025.flux.nc"
data1 = xr.open_dataset(file1)

lat1 = data1['lat']
lon1 = data1['lon']
time1 = data1['mtime']
pco2_1  = data1['co2flux_ocean']
area = data1['dxyp']

co2_flux = pco2_1/area

co2_flux = co2_flux.where(co2_flux!= 0) 

co2_flux =co2_flux*10**15
#%% Area selection
eio= co2_flux.sel(lat = slice(5,22.5),lon = slice(45,65),mtime = slice('1960-01-01','2024-12-31'))

#%% eio

years = np.arange(1960,2025,1)

#%%
eio_1 = eio.resample(mtime = 'QS-DEC').mean(dim=('mtime'))
eio_1 = eio_1.isel(mtime=slice(1, None)) 
eio_d = eio_1.isel(mtime = eio_1.mtime.dt.month == 12)
eio_dclim = eio_d.mean(dim =('lat','lon'))
#%% EIO
fig, axs = plt.subplots(2, 2, figsize=(25,17),dpi = 300)

ax1 = axs[0,0]
x =eio_dclim.mtime.dt.year.values
y = eio_dclim.values
coeff_djf = np.polyfit(x,y,1)
trend_djf = np.polyval(coeff_djf,years)
trend_djf_stat = stats.linregress(x,y)
p_val_djf = trend_djf_stat.pvalue
trend_label_djf = f"Trend ({coeff_djf[0]:.2e} gC/$m^2$yr, p={p_val_djf:.3f})"

ax1.plot(years, eio_dclim, marker='o', color='b', label='Spatially Averaged Mean')
ax1.plot(years,trend_djf, color='red', linestyle='--',linewidth=3, label=trend_label_djf)
ax1.set_title('Mean pCO$_2$ Trend (DJF)', fontsize=24)
ax1.tick_params(axis='both', which='major', labelsize=24)
ax1.grid(True, linestyle=':', alpha=0.7)
ax1.tick_params(axis='both', which='major', labelsize=24)
ax1.legend(fontsize = 20)



other_seasons = [[3, 4, 5], [6, 7, 8, 9], [10, 11]]
season_names = ['MAM', 'JJAS', 'ON']
axes_to_plot = [axs[0, 1], axs[1, 0], axs[1, 1]] 

for i in range(len(other_seasons)):
    ax = axes_to_plot[i] 
    months = other_seasons[i]
    name = season_names[i]

    # B. Use your seasonal calculation logic
    eio_season = eio.sel(mtime=eio.mtime.dt.month.isin(months))
    eio_season = eio_season.resample(mtime='Y').mean(dim='mtime')
    eio_season = eio_season.dropna(dim='mtime', how='all') # Clean empty years
    eio_clim = eio_season.mean(dim=('lat', 'lon'))
    
    # C. Get years and data
    years = eio_clim.mtime.dt.year.values
    data = eio_clim.values

    # D. Calculate trend
    coeff_eio = np.polyfit(years, data, 1)
    trend_line_eio = np.polyval(coeff_eio, years)
    trend_stats = stats.linregress(years, data)
    p_val = trend_stats.pvalue

# Update Label: Include p-value formatted to 3 decimal places
    trend_label = f"Trend ({coeff_eio[0]:.2e} gC/$m^2$yr, p={p_val:.3f})"
    # E. Plot on the correct axis 'ax' (NOT 'plt.')
    ax.plot(years, data, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')
    ax.plot(years, trend_line_eio, color='red', linestyle='--', linewidth=4, label=trend_label)
    
    ax.set_title(f'Mean CO$_2$ flux Trend ({name})', fontsize=25)  
    ax.grid(True, linestyle=':', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.legend(fontsize = 20)

# --- 5. Final Layout ---
axs[0, 0].set_ylim(-4,13)
axs[0, 1].set_ylim(-4,13)
axs[1, 1].set_ylim(-4,13)
axs[1,0].set_ylim(42,60)


# Y-Labels for Plot 1 (Top-Left) and Plot 3 (Bottom-Left)
axs[0, 0].set_ylabel('gC/$m^2$yr', fontsize=24)
axs[1, 0].set_ylabel('gC/$m^2$yr', fontsize=24)

# X-Labels for Plot 3 (Bottom-Left) and Plot 4 (Bottom-Right)
axs[1, 0].set_xlabel('Years', fontsize=24)
axs[1, 1].set_xlabel('Years', fontsize=24)
plt.suptitle('Mean sea- air CO$_2$ flux Time Series for the Northwestern Indian Ocean (5°N-22.5°N, 45°E-65°E)', fontsize=26)
plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust rect to make room for suptitle
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/co2_flux/trend/nwio_CO2F',dpi=600)
plt.show()

#%%
