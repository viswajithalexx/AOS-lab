#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 22:44:01 2025

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
#%% Area selection
eio = pco2_1.sel(lat = slice(-5,5),lon = slice(30,96))
so = pco2_1.sel(lat = slice(-0,24),lon = slice(30,70))
nas = pco2_1.sel(lat = slice(20,30),lon = slice(55,72))
su = pco2_1.sel(lat = slice(-5,10),lon = slice(92,110))
#%% nas_djf_variable
'DJF through resample will calculate DJF of 1957 by (1957 dec+1957 jan+ 1957 Feb) it should be 1957 Dec+1958 Jan and Feb'

nas_q = nas.resample(mtime='QS-DEC').mean(dim=('mtime')) #average of 3 months,for 1957 to 2024
nas_q = nas_q.isel(mtime=slice(1, None)) 
nas_djf = nas_q.sel(mtime=nas_q.mtime.dt.month == 12) # only taking the djf bin
nas_djf_clim = nas_djf.mean(dim=('lat','lon'))
#%% nas_MAM_variable
nas_mam = nas.sel(mtime=nas.mtime.dt.month.isin([3, 4, 5])) #same year no issue using resample
nas_mam = nas_mam.resample(mtime='Y').mean(dim='mtime')
nas_mam_clim = nas_mam.mean(dim = ('lat','lon'))
#%% so-jjas_variable
so_jjas = nas.sel(mtime=so.mtime.dt.month.isin([6,7,8,9]))
so_jjas = so_jjas.resample(mtime='Y').mean(dim='mtime')
so_jjas_clim = so_jjas.mean(dim = ('lat','lon'))
#%% 
su_on = su.sel(mtime=so.mtime.dt.month.isin([10,11]))
su_on = su_on.resample(mtime='Y').mean(dim='mtime')
su_on_clim = su_on.mean(dim = ('lat','lon'))
#%% eio_variable
eio_all = eio.resample(mtime = 'Y').mean(dim ='mtime')
eio_clim = eio_all.mean(dim = ('lat','lon'))
years = np.arange(1957,2025,1)
#%% trend calculation
from scipy import stats
import numpy as np

y_data = su_on_clim.values
x_data = su_on_clim.mtime.dt.year.values
trend_stats = stats.linregress(x_data, y_data)

# 4. Print the results
print(f"Slope (Trend): {trend_stats.slope:.2f} units per year")
print(f"R-value: {trend_stats.rvalue:.2f}")
print(f"P-value: {trend_stats.pvalue:.4f}")
#%% plot_eio_clim
import matplotlib.pyplot as plt
import numpy as np


coeffs_eio = np.polyfit(years, eio_clim, 1)

trend_line_eio = np.polyval(coeffs_eio, years)

plt.figure(figsize=(12, 6), dpi=200)

plt.scatter(years, eio_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')

plt.plot(years, trend_line_eio, color='red', linestyle='--', linewidth=2, label='Linear Trend 1.51 $\mu$atm/year ')

# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('pCO$_2$ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.ylim(300, 470)
plt.title('Mean pCO$_2$ Time Series for the Equatorial Indian Ocean (5°S-5°N, 30°E-96°E)', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend() # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
plt.tight_layout()
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_EIO.png', dpi=500, bbox_inches='tight') 
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_EIO.tiff', dpi=500, bbox_inches='tight')
plt.show()
#%% plot_so_jjas_clim

coeffs_so = np.polyfit(years, so_jjas_clim, 1)

trend_so = np.polyval(coeffs_so, years)

# --- 2. Plot the original data and the trend line ---
plt.figure(figsize=(12, 6), dpi=200)

# Plot your original data
plt.scatter(years, so_jjas_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')

# ➡️ ADDED: Plot the new trend line ⬅️
plt.plot(years, trend_so, color='red', linestyle='--', linewidth=2, label='Linear Trend 1.27 $\mu$atm/year' )

# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('pCO$_2$ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.ylim(300, 470)
plt.title('Seasonal (JJAS) pCO$_2$ Time Series for the Northwestern Indian Ocean (0°-24°N, 30°-70°E)', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend() # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
plt.tight_layout()
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_SO_jjas.png', dpi=500, bbox_inches='tight') 
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_SO_jjas.tiff', dpi=500, bbox_inches='tight') 

plt.show()
#%% plot_nas_mam_clim
coeffs_nas_mam = np.polyfit(years,nas_mam_clim, 1)



trend_nas_mam = np.polyval(coeffs_nas_mam, years)


plt.figure(figsize=(12, 6), dpi=200)


plt.scatter(years, nas_mam_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')

# ➡️ ADDED: Plot the new trend line ⬅️
plt.plot(years, trend_nas_mam, color='red', linestyle='--', linewidth=2, label='Linear Trend 1.31 $\mu$atm/year')

# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('Spatially Averaged pCO$_2$ ($\mu$atm)', fontsize=12) 
plt.ylim(300, 470)
plt.title('Seasonal (MAM) pCO$_2$ Time Series for the Northern Arabian Sea (20°N-30°N, 55°E-72°E)', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend() 
plt.tight_layout()
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_NAS_mam.png', dpi=500, bbox_inches='tight') 
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_NAS_mam.tiff', dpi=500, bbox_inches='tight') 

plt.show()
#%% plot nas_djf
coeffs_nas_djf = np.polyfit(years,nas_djf_clim, 1)

trend_nas_djf = np.polyval(coeffs_nas_djf, years)

plt.figure(figsize=(12, 6), dpi=200)

plt.scatter(years, nas_djf_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')

plt.plot(years, trend_nas_djf, color='red', linestyle='--', linewidth=2, label='Linear Trend 1.28 $\mu$atm/year')

# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('pCO$_2$ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.ylim(300, 470)
plt.title('Seasonal (DJF) pCO$_2$ Time Series for the Northern Arabian Sea (20°N-30°N, 55°E-72°E)', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend() # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
plt.tight_layout()
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_NAS_DJF.png', dpi=500, bbox_inches='tight') 
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_NAS_DJF.tiff', dpi=500, bbox_inches='tight') 

plt.show()
#%%
coeffs_su = np.polyfit(years, su_on_clim, 1)

trend_su = np.polyval(coeffs_su, years)-
plt.figure(figsize=(12, 6), dpi=200)
plt.scatter(years, su_on_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')
plt.plot(years, trend_su, color='red', linestyle='--', linewidth=2, label='Linear Trend 1.21 $\mu$atm/year' )
plt.xlabel('Year', fontsize=12)
plt.ylabel('pCO$_2$ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.ylim(300, 470)
plt.title('Seasonal (ON) pCO$_2$ Time Series for the Eastern Indian Ocean (5°S-10°N, 92°E-110°E)', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend() # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
plt.tight_layout()
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_su_on.png', dpi=500, bbox_inches='tight') 
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/trend_pco2_su_on.tiff', dpi=500, bbox_inches='tight') 

plt.show()
#%