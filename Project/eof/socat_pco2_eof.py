#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 21:19:34 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
from eofs.xarray import Eof
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean.cm as cm

file1 = "/home/bobco-08/24cl05012/CO2/data/oc_v2025.pCO2.nc"
data1 = xr.open_dataset(file1)

lat1 = data1['lat']
lon1 = data1['lon']
time1 = data1['mtime']
pco2_1  = data1['pCO2']

pco2_1 = pco2_1.where(pco2_1 != 0) 
#%%
socat = pco2_1.sel(lat = slice(-30,30),lon = slice(30,110))
#%%
da = socat
#%%
da_mon = da.groupby('mtime.month')
da_clim = da_mon.mean()
da_rm = da_mon - da_clim
#%%




#%%
# --- Extract the first 3 EOFs and PCs ---
eofs1  = solver1.eofs(neofs=3)   # spatial patterns
pcs1   = solver1.pcs(npcs=3, pcscaling=1)     # time series
frac1  = solver1.varianceFraction(neigs=3)    # variance explained
# --- Print explained variance ---
for i in range(3):
    print(f"EOF{i+1} explains {frac1[i].values*100:.2f}% variance")

#%%
    plt.figure(figsize=(15,10), dpi=200) 
    plt.plot(pcs1.time, pcs1.isel(mode=3)) 
    plt.xlabel('YEAR')
    plt.ylabel('\u00b5atm')
    plt.title(f'SOCAT pCO₂ EOF Mode {i+1} ({frac1[i].values*100:.2f}% variance explained)')
    plt.grid(True)
    plt.show()
#%%
da_mon = da.groupby('mtime.month')
da_clim = da_mon.mean()
da_rm = da_mon - da_clim

#%% detrend

# Fit a linear trend along 'mtime'
trend = da_rm.polyfit(dim='mtime', deg=1)  # degree=1 → linear
#%%

# Evaluate trend at each time step
da_trend = xr.polyval(da_rm['mtime'], trend.polyfit_coefficients)
#%%
# Subtract trend to get detrended data
da_detrended = da_rm - da_trend
#%%
da_ocean = da_rm.where(~da_rm.isnull().all(dim = ('mtime')))
#%%
da_ocean_mo = da_ocean.resample(mtime ='M').mean()

da_ocean_mo = da_ocean_mo.rename({'mtime':'time'})
solver = Eof(da_ocean_mo)

# --- Extract the first 3 EOFs and PCs ---
eofs  = solver.eofs(neofs=3)   # spatial patterns
pcs   = solver.pcs(npcs=3, pcscaling=1)     # time series
frac  = solver.varianceFraction(neigs=3)    # variance explained
# --- Print explained variance ---
for i in range(3):
    print(f"EOF{i+1} explains {frac[i].values*100:.2f}% variance")
#%%
for i in range(3):
    plt.figure(figsize = (15,10),dpi = 300)
    ax = plt.axes(projection = ccrs.PlateCarree())
    ax.coastlines()
    levels = np.linspace(-0.18,0.18,40)
    im =  plt.contourf(eofs.lon,eofs.lat,eofs[i],levels,cmap = cm.balance)
    ax.add_feature(cf.LAND,color = 'grey',zorder =10)
    ax.gridlines(visible=True,draw_labels=True)
    cbar = plt.colorbar(im)
    cbar.set_label('\u00b5atm')
    ax.set_title(f'SOCAT pCO₂ EOF Mode {i+1} ({frac[i].values*100:.2f}% variance explained)', pad = 30)
#%% 
for i in range(3):
    plt.figure(figsize=(15,10), dpi=200)
    plt.plot(pcs.time, pcs.isel(mode=i)) 
    plt.ylim(-400,400)
    plt.xlabel('YEAR')
    plt.ylabel('\u00b5atm')
    plt.title(f'SOCAT pCO₂ EOF Mode {i+1} ({frac[i].values*100:.2f}% variance explained)')
    plt.grid(True)
    plt.show()

#%%

da_ocean_yr = da_ocean.groupby('mtime.year').mean()

da_ocean_yr = da_ocean_yr.rename({'year':'time'})
solver = Eof(da_ocean_yr)

# --- Extract the first 3 EOFs and PCs ---
eofs_y  = solver.eofs(neofs=3)   # spatial patterns
pcs_y   = solver.pcs(npcs=3, pcscaling=0)     # time series
frac_y  = solver.varianceFraction(neigs=3)    # variance explained

# --- Print explained variance ---
for i in range(3):
    print(f"EOF{i+1} explains {frac[i].values*100:.2f}% variance")
#%%
for i in range(3):
    plt.figure(figsize = (15,10),dpi = 300)
    ax = plt.axes(projection = ccrs.PlateCarree())
    ax.coastlines()
    
    im =  plt.contourf(eofs.lon,eofs.lat,eofs_y[i],levels,cmap = cm.balance)
    ax.add_feature(cf.LAND,color = 'grey',zorder =10)
    ax.gridlines(visible=True,draw_labels=True)
    cbar = plt.colorbar(im)
    cbar.set_label('\u00b5atm')
    ax.set_title(f'SOCAT pCO₂ EOF Mode {i+1} ({frac[i].values*100:.2f}% variance explained)', pad = 30)
#%% 
for i in range(3):
    plt.figure(figsize=(15,10), dpi=200)
    plt.plot(pcs_y.time, pcs_y.isel(mode=i)) 
    plt.ylim(-400,400)
    plt.xlabel('YEAR')
    plt.ylabel('Standardized Units')
    plt.title(f'SOCAT pCO₂ EOF Mode {i+1} ({frac[i].values*100:.2f}% variance explained)')
    plt.grid(True)
    plt.show()

#%%
