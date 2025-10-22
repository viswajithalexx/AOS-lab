#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 19:51:42 2025

@author: bobco-08
"""


import xarray as xr
import numpy as np
from eofs.xarray import Eof
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean.cm as cm

file0 = '/home/bobco-08/24cl05012/CO2/data/pco2_climatology_taka.nc'
file1 = "/home/bobco-08/24cl05012/CO2/data/oc_v2024E.pCO2.nc"
file2 = '/home/bobco-08/24cl05012/CO2/data/pCO2-Corrected_INCOIS-BIO-ROMS_v2.nc'

data0 = xr.open_dataset(file0)
data1 = xr.open_dataset(file1)
data2 = xr.open_dataset(file2)

lat0 = data0['lat']
lon0 = data0['lon']
mon = data0['month']
pco2_0 = data0['pco2_sw']

lat1 = data1['lat']
lon1 = data1['lon']
time1 = data1['mtime']
pco2_1  = data1['pCO2']

pco2_1 = pco2_1.where(pco2_1 != 0) 

lat2 = data2['LAT']
lon2 = data2['LON']
time2 = data2['TIME']
pco2_2 = data2['pCO2_Original']
#%%

clim = pco2_0.sel(lat = slice(0,30),lon = slice(38,110))

incois = pco2_2.sel(LAT = slice(0,30),LON = slice(38,110),TIME = slice('1980-01-01','2007-12-31'))

socat = pco2_1.sel(lat = slice(0,30),lon = slice(38,100),mtime = slice('1980-01-01','2007-12-31'))

#%% 3 month running mean

da = socat
#running mean
da_monthly = da.resample(mtime='M').mean()
da_rm = da_monthly.rolling(mtime=3, center=True).mean()

#spatial plot of running mean
plt.figure(figsize = (15,10),dpi = 200)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
levels = np.linspace(250,500,60)
im =  plt.contourf(da_rm.lon,da_rm.lat,da_rm.mean(dim=('mtime')),levels,cmap = cm.balance)
ax.add_feature(cf.LAND,color = 'grey',zorder =10)
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
cbar.set_label('\u00b5atm')
ax.set_title('3-months running mean of pCO₂ IO ', pad = 30)
#%% detrend

# Fit a linear trend along 'mtime'
trend = da_rm.polyfit(dim='mtime', deg=1)  # degree=1 → linear

# Evaluate trend at each time step
da_trend = xr.polyval(da_rm['mtime'], trend.polyfit_coefficients)

# Subtract trend to get detrended data
da_detrended = da_rm - da_trend

da_detrended = da_detrended.interpolate_na(dim= 'mtime',method = 'linear', fill_value="extrapolate")

da_ocean = da_detrended.where(~da_detrended.isnull().all(dim = ('mtime')))


solver = Eof(da_ocean)

# --- Extract the first 3 EOFs and PCs ---
eofs  = solver.eofsAsCorrelation(neofs=3)   # spatial patterns
pcs   = solver.pcs(npcs=3, pcscaling=1)     # time series
frac  = solver.varianceFraction(neigs=3)    # variance explained

# --- Print explained variance ---
for i in range(3):
    print(f"EOF{i+1} explains {frac[i].values*100:.2f}% variance")
#%%

plt.figure(figsize = (15,10),dpi = 200)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
im =  plt.contourf(eofs.lon,eofs.lat,eofs[0],levels = 20,cmap = 'RdBu_r')
ax.add_feature(cf.LAND,color = 'grey',zorder =10)
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
# cbar.set_label('\u00b5atm')
ax.set_title('EOF mode 1 ', pad = 30)

#%%

plt.figure(figsize = (15,10),dpi = 200)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
im =  plt.contourf(eofs.lon,eofs.lat,eofs[1],levels = 20,cmap = 'RdBu_r')
ax.add_feature(cf.LAND,color = 'grey',zorder =10)
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
# cbar.set_label('\u00b5atm')
ax.set_title('EOF mode 2 ', pad = 30)
#%%

plt.figure(figsize = (15,10),dpi = 200)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
im =  plt.contourf(eofs.lon,eofs.lat,eofs[2],levels = 20,cmap = 'RdBu_r')
ax.add_feature(cf.LAND,color = 'grey',zorder =10)
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
# cbar.set_label('\u00b5atm')
ax.set_title('EOF mode 3 ', pad = 30)