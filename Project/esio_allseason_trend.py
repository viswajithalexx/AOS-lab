#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  8 02:04:32 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


file1 = "/home/bobco-08/24cl05012/CO2/data/oc_v2025.pCO2.nc"
data1 = xr.open_dataset(file1)

lat1 = data1['lat']
lon1 = data1['lon']
time1 = data1['mtime']
pco2_1  = data1['pCO2']

pco2_1 = pco2_1.where(pco2_1 != 0) 
#%% Area selection
esio = pco2_1.sel(lat = slice(-5,10),lon = slice(92,110))
#%% eio
years = np.arange(1957,2025,1)

#%%
esio_1 = esio.resample(mtime = 'QS-DEC').mean(dim=('mtime'))
esio_1 = esio_1.isel(mtime=slice(1, None)) 
esio_d = esio_1.isel(mtime = esio_1.mtime.dt.month == 12)
esio_dclim = esio_d.mean(dim =('lat','lon'))
#%% EIO
fig, axs = plt.subplots(2, 2, figsize=(20,18),dpi = 600)

ax1 = axs[0,0]
x = esio_dclim.mtime.dt.year.values
y = esio_dclim.values
coeff_djf = np.polyfit(x,y,1)
trend_djf = np.polyval(coeff_djf,years)
trend_djf_stat = stats.linregress(x,y)
trend_label_djf = f"Trend ({coeff_djf[0]:.2f} $\mu$atm/yr)"

ax1.plot(years, esio_dclim, marker='o', color='b', label='Spatially Averaged Mean')
ax1.plot(years,trend_djf, color='red', linestyle='--', linewidth=3, label=trend_label_djf)
ax1.set_ylim(280,480)
ax1.set_xlabel('Years', fontsize=24)
ax1.set_ylabel('$\mu$atm', fontsize=24)
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
    esio_season = esio.sel(mtime=esio.mtime.dt.month.isin(months))
    esio_season = esio_season.resample(mtime='Y').mean(dim='mtime')
    esio_season = esio_season.dropna(dim='mtime', how='all') # Clean empty years
    esio_clim = esio_season.mean(dim=('lat', 'lon'))
    
    # C. Get years and data
    years = esio_clim.mtime.dt.year.values
    data = esio_clim.values

    # D. Calculate trend
    coeff_eio = np.polyfit(years, data, 1)
    trend_line_eio = np.polyval(coeff_eio, years)
    trend_stats = stats.linregress(years, data)
    trend_label = f"Trend ({coeff_eio[0]:.2f} $\mu$atm/yr)"
    ax.set_ylim(280,480)
    # E. Plot on the correct axis 'ax' (NOT 'plt.')
    ax.plot(years, data, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')
    ax.plot(years, trend_line_eio, color='red', linestyle='--', linewidth=3, label=trend_label)  
    ax.set_xlabel('Years', fontsize = 24)
    ax.set_ylabel('$\mu$atm',fontsize = 24)
    ax.set_title(f'Mean pCO$_2$  Trend ({name})', fontsize=25)
    ax.grid(True, linestyle=':', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.legend(fontsize = 20)

# --- 5. Final Layout ---
plt.suptitle('Mean pCO$_2$ Time Series for the Eastern Indian Ocean (5°S-10°N, 92°E-110°E)', fontsize=26)
plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust rect to make room for suptitle
plt.show()

