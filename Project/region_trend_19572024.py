#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 23:00:35 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
from eofs.xarray import Eof
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean.cm as cm
import cmaps

file1 = "/home/bobco-08/24cl05012/CO2/data/oc_v2025.pCO2.nc"
data1 = xr.open_dataset(file1)

lat1 = data1['lat']
lon1 = data1['lon']
time1 = data1['mtime']
pco2_1  = data1['pCO2']

pco2_1 = pco2_1.where(pco2_1 != 0) 
#%%  region selection

eio = pco2_1.sel(lat = slice(-6.5,5),lon = slice(49,92))
was = pco2_1.sel(lat = slice(5,22.5),lon = slice(45,65))
nas = pco2_1.sel(lat = slice(22.5,28),lon = slice(56,70))
esio = pco2_1.sel(lat = slice(-6.6,8),lon = slice(92,109))
#%% time series plot 1957 - 2024

eio_m = eio.groupby('mtime.month').mean()
eio_t = eio_m.mean(dim = ('lat','lon'))
was_m  = was.groupby('mtime.month').mean()
was_t = was_m.mean(dim=('lat','lon')) 
nas_m  = nas.groupby('mtime.month').mean()
nas_t = nas_m.mean(dim=('lat','lon'))
esio_m  = esio.groupby('mtime.month').mean()
esio_t  = esio_m.mean(dim=('lat','lon'))

#%% region wise plot

months = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
tseries = [eio_t ,was_t,nas_t,esio_t]
region =['the Equatorial Indian Ocean (5°S-5°N, 30°E-96°E)','the Northwestern Indian Ocean (0°-24°N, 30°-70°E)','the Northern Arabian Sea (20°N-30°N, 55°E-72°E)','the Eastern Indian Ocean (5°S-10°N, 92°E-110°E)']

for j in range(len(tseries)):
        plt.figure(figsize=(12, 6), dpi=200)

        plt.plot(months,tseries[j], marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')
        # Add labels and title
        plt.xlabel('Months', fontsize=12)
        plt.ylabel('pCO$_2$ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
        plt.title(f'Monthly Climatology of pCO$_2$ for {region[j]}', fontsize=14)
        plt.grid(True, linestyle=':', alpha=0.7)
        plt.legend() # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
        plt.tight_layout()
        plt.show()
#%%   common plot      

plt.figure(figsize=(12, 6), dpi=600)

plt.plot(months,eio_t, marker='o', linestyle='-', color='r', label='EIO (6.5°S-5°N, 49°E-92°E)')
plt.plot(months,was_t, marker='o', linestyle='-', color='b', label='NWIO (5°N-22.5°N, 45°E-65°E)')
plt.plot(months,nas_t, marker='o', linestyle='-', color='g', label='NAS (22.5°N-28°N, 56°E-70°E)')
plt.plot(months,esio_t, marker='o', linestyle='-', color='purple', label='ESIO (6.6°S-8°N, 92°E-109°E)')
# Add labels and title
plt.ylim(310,460)
plt.xlabel('Months', fontsize=20)
plt.ylabel('pCO$_2$ ($\mu$atm)', fontsize=18) # Assuming pCO2 units
plt.title('Monthly Climatology of pCO$_2$ of various Indian Ocean regions', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(fontsize =12)
plt.tick_params(axis='both', which='major', labelsize=14) # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
plt.tight_layout()
plt.show()
#%% 