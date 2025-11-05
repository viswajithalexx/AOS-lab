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
#%%

eio = pco2_1.sel(lat = slice(-5,5),lon = slice(30,96))
so = pco2_1.sel(lat = slice(-0,24),lon = slice(30,70))
nas = pco2_1.sel(lat = slice(20,30),lon = slice(55,72))

#%%

plt.figure(figsize =(15,10),dpi = 300)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines(zorder = 15)
levels = np.linspace(300,650,25)

im =  plt.contourf(nas.lon, nas.lat, nas.mean(dim =('mtime')),levels,cmap = cmaps.cmp_b2r)
my_contour_levels = np.arange(300,490, 20)

ax.add_feature(cf.LAND,color = 'grey',zorder =12)
ax.gridlines(visible= False ,draw_labels=True)
cbar_ticks = np.arange(370, 650, 20)
cbar = plt.colorbar(im,ticks=cbar_ticks)
cbar.set_label('\u00b5atm')
plt.show()
#%%

nas_djf = nas.sel(mtime=nas.mtime.dt.month.isin([12, 1, 2]))

nas_djf = nas_djf.resample(mtime='Y').mean(dim='mtime')

nas_djf_clim = nas_djf.mean(dim = ('lat','lon'))
#%%

nas_mam = nas.sel(mtime=nas.mtime.dt.month.isin([3, 4, 5]))

nas_mam = nas_mam.resample(mtime='Y').mean(dim='mtime')

nas_mam_clim = nas_mam.mean(dim = ('lat','lon'))
#%%
so_jjas = nas.sel(mtime=so.mtime.dt.month.isin([6,7,8,9]))

so_jjas = so_jjas.resample(mtime='Y').mean(dim='mtime')

so_jjas_clim = so_jjas.mean(dim = ('lat','lon'))
#%%
eio_all = eio.resample(mtime = 'Y').mean(dim ='mtime')

eio_clim = eio_all.mean(dim = ('lat','lon'))
#%% nas_djf
# Create the line plot
years = np.arange(1957,2025,1)
plt.figure(figsize=(12, 6), dpi=200)
plt.plot(years, nas_djf_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')
plt.ylim(310,450)
# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('Spatially Averaged pCO₂ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.title('Spatially Averaged DJF Mean pCO₂ Time Series', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.tight_layout()
#%% nas_mam
plt.figure(figsize=(12, 6), dpi=200)
plt.plot(years, nas_mam_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')
plt.ylim(310,450)
# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('Spatially Averaged pCO₂ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.title('Spatially Averaged MAM Mean pCO₂ Time Series', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.tight_layout()
#%%so_jjas
plt.figure(figsize=(12, 6), dpi=200)
plt.plot(years, so_jjas_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')
plt.ylim(310,450)
# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('Spatially Averaged pCO₂ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.title('Spatially Averaged MAM Mean pCO₂ Time Series', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.tight_layout()
#%% eio
plt.figure(figsize=(12, 6), dpi=200)
plt.plot(years, eio_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')
plt.ylim(310,450)
# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('Spatially Averaged pCO₂ ($\mu$atm)', fontsize=12) # Assuming pCO2 units

plt.title('box average of EIO', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.tight_layout()
#%%
from scipy import stats
import numpy as np

# --- Assuming 'nas_djf_clim' is your 1D DataArray (time) ---
# nas_djf_clim = nas_djf.mean(dim=('lat','lon'))

# 1. Get the Y-values (the data)
y_data = nas_djf_clim.values

# 2. Get the X-values (the time component)
#    It's best to use the years as numbers
x_data = nas_djf_clim.mtime.dt.year.values

# 3. Calculate the linear regression statistics
#    This returns slope, intercept, r-value, p-value, and std-error
trend_stats = stats.linregress(x_data, y_data)

# 4. Print the results
print(f"Slope (Trend): {trend_stats.slope:.2f} units per year")
print(f"R-value: {trend_stats.rvalue:.2f}")
print(f"P-value: {trend_stats.pvalue:.4f}")
#%%
import matplotlib.pyplot as plt
import numpy as np

# --- Assumed variables ---
# You must have 'years' (1D array) and 'eio_clim' (1D array) defined.
# For example:
# years = np.arange(1980, 2020)
# eio_clim = 350 + np.arange(40) * 0.5 + np.random.randn(40) * 5
# -------------------------

# --- 1. Calculate the trend line ---
# Use np.polyfit() to find the slope and intercept of a 1st-degree (linear) fit
# coeffs[0] will be the slope, coeffs[1] will be the intercept
coeffs = np.polyfit(years, eio_clim, 1)

# Use np.polyval() to create the y-values for the trend line
# This evaluates the linear equation at all 'years' points
trend_line = np.polyval(coeffs, years)

# --- 2. Plot the original data and the trend line ---
plt.figure(figsize=(12, 6), dpi=200)

# Plot your original data
plt.plot(years, eio_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')

# ➡️ ADDED: Plot the new trend line ⬅️
plt.plot(years, trend_line, color='red', linestyle='--', linewidth=2, label='Linear Trend')

# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('Spatially Averaged pCO₂ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.ylim(320, 470)
plt.title('box average of EIO', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend() # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
plt.tight_layout()
plt.show()
#%%

coeffs = np.polyfit(years, so_jjas_clim, 1)

# Use np.polyval() to create the y-values for the trend line
# This evaluates the linear equation at all 'years' points
trend_so_jjas_clim = np.polyval(coeffs, years)

# --- 2. Plot the original data and the trend line ---
plt.figure(figsize=(12, 6), dpi=200)

# Plot your original data
plt.plot(years, so_jjas_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')

# ➡️ ADDED: Plot the new trend line ⬅️
plt.plot(years, trend_so_jjas_clim, color='red', linestyle='--', linewidth=2, label='Linear Trend')

# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('Spatially Averaged pCO₂ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.ylim(320, 470)
plt.title('Spatially Averaged MAM Mean pCO₂ Time Series', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend() # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
plt.tight_layout()
plt.show()
#%%
coeffs = np.polyfit(years,nas_mam_clim, 1)

# Use np.polyval() to create the y-values for the trend line
# This evaluates the linear equation at all 'years' points
trend_nas_mam_clim = np.polyval(coeffs, years)

# --- 2. Plot the original data and the trend line ---
plt.figure(figsize=(12, 6), dpi=200)

# Plot your original data
plt.plot(years, nas_mam_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')

# ➡️ ADDED: Plot the new trend line ⬅️
plt.plot(years, trend_nas_mam_clim, color='red', linestyle='--', linewidth=2, label='Linear Trend')

# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('Spatially Averaged pCO₂ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.ylim(320, 470)
plt.title('Spatially Averaged MAM Mean pCO₂ Time Series', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend() # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
plt.tight_layout()
plt.show()
#%%
coeffs = np.polyfit(years,nas_djf_clim, 1)

# Use np.polyval() to create the y-values for the trend line
# This evaluates the linear equation at all 'years' points
trend_nas_djf_clim = np.polyval(coeffs, years)

# --- 2. Plot the original data and the trend line ---
plt.figure(figsize=(12, 6), dpi=200)

# Plot your original data
plt.plot(years, nas_djf_clim, marker='o', linestyle='-', color='b', label='Spatially Averaged Mean')

# ➡️ ADDED: Plot the new trend line ⬅️
plt.plot(years, trend_nas_djf_clim, color='red', linestyle='--', linewidth=2, label='Linear Trend')

# Add labels and title
plt.xlabel('Year', fontsize=12)
plt.ylabel('Spatially Averaged pCO₂ ($\mu$atm)', fontsize=12) # Assuming pCO2 units
plt.ylim(320, 470)
plt.title('Spatially Averaged DJF Mean pCO₂ Time Series', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend() # This will now show both 'Spatially Averaged Mean' and 'Linear Trend'
plt.tight_layout()
plt.show()