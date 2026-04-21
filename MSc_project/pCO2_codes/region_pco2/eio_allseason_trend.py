#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  8 00:10:46 2025

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmaps


file = "/home/bobco-08/24cl05012/CO2/data/data_1/cmems/bio_var/cmems_obs-mob_glo_bgc-car_my_irr-i_1775027250539_30E-119.88E_29.88S-29.88N.nc"

data = xr.open_dataset(file)

lat =  data['latitude']
lon = data['longitude']
time = data['time']
pco2 = data['spco2']
# surface douwnward flux of total carbon
pco2v = pco2.values
# for IO

pco2_io = pco2.sel(time = slice('1994-01-01','2024-12-01'))
#%% djf mean pco2

pco2_djf = pco2_io.resample(time = 'QS-DEC').mean(dim=('time'))
pco2_djf = pco2_djf.isel(time=slice(1, None))
pco2_djf = pco2_djf.isel(time = pco2_djf.time.dt.month == 12)

djf_mean = pco2_djf.mean(dim ='time')

#%% mam mean pco2

pco2_mam = pco2_io.sel(time= pco2_io.time.dt.month.isin([3,4,5]))
pco2_mam = pco2_mam.resample(time='YE').mean(dim='time')
mam_mean = pco2_mam.mean(dim='time')

#%% jjas mean pco2
pco2_jjas = pco2_io.sel(time= pco2_io.time.dt.month.isin([6,7,8,9]))
pco2_jjas = pco2_jjas.resample(time='YE').mean(dim='time')
jjas_mean = pco2_jjas.mean(dim='time')

#%% on mean pco2

pco2_on = pco2_io.sel(time= pco2_io.time.dt.month.isin([10,11]))
pco2_on = pco2_on.resample(time='YE').mean(dim='time')
on_mean = pco2_on.mean(dim='time')

#%% EIO
fig, axs = plt.subplots(2, 2, figsize=(20,18),dpi = 600)

ax1 = axs[0,0]
x =eio_dclim.mtime.dt.year.values
y = eio_dclim.values
coeff_djf = np.polyfit(x,y,1)
trend_djf = np.polyval(coeff_djf,years)
trend_djf_stat = stats.linregress(x,y)
trend_label_djf = f"Trend ({coeff_djf[0]:.2f} $\mu$atm/yr)"

ax1.plot(years, eio_dclim, marker='o', color='b', label='Spatially Averaged Mean')
ax1.plot(years,trend_djf, color='red', linestyle='--',linewidth=3, label=trend_label_djf)
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
    trend_label = f"Trend ({coeff_eio[0]:.2f} $\mu$atm/yr)"
       
    # E. Plot on the correct axis 'ax' (NOT 'plt.')
    ax.plot(years, data, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')
    ax.plot(years, trend_line_eio, color='red', linestyle='--', linewidth=4, label=trend_label)
    ax.set_ylim(280,480)
    ax.set_xlabel('Years', fontsize = 24)
    ax.set_ylabel('$\mu$atm',fontsize = 24)
    ax.set_title(f'Mean pCO$_2$  Trend ({name})', fontsize=25)  
    ax.grid(True, linestyle=':', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.legend(fontsize = 20)

# --- 5. Final Layout ---
plt.suptitle('Mean pCO$_2$ Time Series for the Equatorial Indian Ocean (6.5°S-5°N, 49°E-92°E)', fontsize=26)
plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust rect to make room for suptitle
plt.show()

