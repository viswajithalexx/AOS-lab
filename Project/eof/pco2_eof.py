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
#%%
data0 = xr.open_dataset(file0)
data1 = xr.open_dataset(file1)
data2 = xr.open_dataset(file2, chunks={'TIME': 12, 'LAT': 50, 'LON': 50})
#%%
lat0 = data0['lat']
lon0 = data0['lon']
mon = data0['month']
pco2_0 = data0['pco2_sw']
#%%
lat1 = data1['lat']
lon1 = data1['lon']
time1 = data1['mtime']
pco2_1  = data1['pCO2']
#%%
pco2_1 = pco2_1.where(pco2_1 != 0) 
#%%
lat2 = data2['LAT']
lon2 = data2['LON']
time2 = data2['TIME']
pco2_2 = data2['pCO2_Original']
#%%

clim = pco2_0.sel(lat = slice(-30,30),lon = slice(30,120))
#%%

incois = pco2_2.sel(LAT = slice(-30,30),LON = slice(30,120),TIME = slice('1980-01-01','2019-12-31'))
#%%
socat = pco2_1.sel(lat = slice(-30,30),lon = slice(30,120),mtime = slice('1980-01-01','2019-12-31'))

#%% 3 month running mean

da = socat
#running mean
da_monthly = da.resample(mtime='M').mean()
da_rm = da_monthly.rolling(mtime=3, center=True).mean()

#%%
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
ax.set_title('3-months running mean of pCO₂ IO socat ', pad = 30)
#%% detrend

# Fit a linear trend along 'mtime'
trend = da_monthly.polyfit(dim='mtime', deg=1)  # degree=1 → linear

# Evaluate trend at each time step
da_trend = xr.polyval(da_monthly['mtime'], trend.polyfit_coefficients)

# Subtract trend to get detrended data
da_detrended = da_monthly - da_trend

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
for i in range(3):
    plt.figure(figsize = (15,10),dpi = 300)
    ax = plt.axes(projection = ccrs.PlateCarree())
    ax.coastlines()
    levels = np.linspace(-1,1,30)
    im =  plt.contourf(eofs.lon,eofs.lat,eofs[i],levels,cmap = 'RdBu_r')
    ax.add_feature(cf.LAND,color = 'grey',zorder =10)
    ax.gridlines(visible=True,draw_labels=True)
    cbar = plt.colorbar(im)
    cbar.set_label('\u00b5atm')
    ax.set_title(f'SOCAT pCO₂ EOF Mode {i+1} ({frac[i].values*100:.2f}% variance explained)', pad = 30)
 #%%
 
for i in range(3):
    plt.figure(figsize=(15,10), dpi=200)
    plt.plot(pcs.mtime, pcs.isel(mode=i)) 
    plt.xlabel('YEAR')
    plt.ylabel('\u00b5atm')
    plt.title(f'SOCAT pCO₂ EOF Mode {i+1} ({frac[i].values*100:.2f}% variance explained)')
    plt.grid(True)
    plt.show()

    
#%% 3 month running mean (climatology)

dc = clim

dc_rm = dc.rolling(month=3, center=True).mean()

#spatial plot of running mean
plt.figure(figsize = (15,10),dpi = 200)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
levels = np.linspace(250,500,60)
im =  plt.contourf(dc_rm.lon,dc_rm.lat,dc_rm.mean(dim=('month')),levels,cmap = cm.balance)
ax.add_feature(cf.LAND,color = 'grey',zorder =10)
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
cbar.set_label('\u00b5atm')
ax.set_title('3-months running mean of pCO₂ IO (climatology)' ,pad = 30)

#%%
# trend = dc_rm.polyfit(dim='month', deg=1)  # degree=1 → linear
# dc_trend = xr.polyval(dc_rm['month'], trend.polyfit_coefficients)
# # Subtract trend to get detrended data
# dc_detrended = dc_rm - dc_trend
# dc_ocean = dc_detrended.where(~dc_detrended.isnull().all(dim = ('month')))

dc_rm = dc_rm.rename({'month':'time'})

dc_rm = dc_rm.interpolate_na(dim= 'time',method = 'linear', fill_value="extrapolate")

solver_c = Eof(dc_rm)

# --- Extract the first 3 EOFs and PCs ---
eofs_c  = solver_c.eofsAsCorrelation(neofs=3)   # spatial patterns
pcs_c   = solver_c.pcs(npcs=3, pcscaling=1)     # time series
frac_c  = solver_c.varianceFraction(neigs=3)    # variance explained

# --- Print explained variance ---
for i in range(3):
    print(f"LDEO v2019 pCO₂ EOF mode {i+1} (explains {frac_c[i].values*100:.2f}% variance)"
)
#%%

for i in range(3):
    plt.figure(figsize = (15,10),dpi = 200)
    ax = plt.axes(projection = ccrs.PlateCarree())
    ax.coastlines()
    im =  plt.contourf(eofs_c.lon,eofs_c.lat,eofs_c[i],levels = 20,cmap = 'RdBu_r')
    ax.add_feature(cf.LAND,color = 'grey',zorder =10)
    ax.gridlines(visible=True,draw_labels=True)
    cbar = plt.colorbar(im)
    # cbar.set_label('\u00b5atm')
    ax.set_title(f'EOF mode {i+1} (explains {frac_c[i].values*100:.2f}% variance) ', pad = 30)
    
#%%
for i in range(3):
    plt.figure(figsize=(15,10), dpi=200)
    plt.plot(pcs_c.isel(mode=i)) 
    plt.xlabel('YEAR')
    plt.ylabel('\u00b5atm')
    plt.title(f'LDEO v2019  pCO₂ EOF Mode {i+1} ({frac[i].values*100:.2f}% variance explained)')
    plt.grid(True)
    plt.show()

#%% 3 month running mean (incois)

di = incois
#%%
di_rm = di.rolling(TIME=3, center=True)
di_rm_m =di_rm.mean(dim=('TIME'))

#%%
#spatial plot of running mean
plt.figure(figsize = (15,10),dpi = 200)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
levels = np.linspace(250,500,60)
im =  plt.contourf(di_rm.LON,di_rm.LAT,di_rm.mean(dim=('TIME')),levels,cmap = cm.balance)
ax.add_feature(cf.LAND,color = 'grey',zorder =10)
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
cbar.set_label('\u00b5atm')
ax.set_title('3-months running mean of pCO₂ IO (INCOIS)' ,pad = 30)
#%%
trend = di_rm.polyfit(dim='TIME', deg=1)  # degree=1 → linear

di_trend = xr.polyval(di_rm['TIME'], trend.polyfit_coefficients)

# Subtract trend to get detrended data
di_detrended = di_rm - di_trend

di_detrended = di_detrended.interpolate_na(dim= 'TIME',method = 'linear', fill_value="extrapolate")

di_ocean = di_detrended.where(~di_detrended.isnull().all(dim = ('TIME')))


solver_i = Eof(di_ocean)

# --- Extract the first 3 EOFs and PCs ---
eofs_i  = solver_i.eofsAsCorrelation(neofs=3)   # spatial patterns
pcs_i   = solver_i.pcs(npcs=3, pcscaling=1)     # time series
frac_i  = solver_i.varianceFraction(neigs=3)    # variance explained

# --- Print explained variance ---
for i in range(3):
    print(f"INCOIS pCO₂ EOF mode {i+1} (explains {frac_i[i].values*100:.2f}% variance)"
)
#%%

for i in range(3):
    plt.figure(figsize = (15,10),dpi = 200)
    ax = plt.axes(projection = ccrs.PlateCarree())
    ax.coastlines()
    im =  plt.contourf(eofs_i.LON,eofs_i.LAT,eofs_i[i],levels = 20,cmap = 'RdBu_r')
    ax.add_feature(cf.LAND,color = 'grey',zorder =10)
    ax.gridlines(visible=True,draw_labels=True)
    cbar = plt.colorbar(im)
    # cbar.set_label('\u00b5atm')
    ax.set_title(f'EOF mode {i+1} (explains {frac_i[i].values*100:.2f}% variance) ', pad = 30)
#%%    

