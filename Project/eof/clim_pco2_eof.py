#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 21:23:18 2025

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

data0 = xr.open_dataset(file0)

lat0 = data0['lat']
lon0 = data0['lon']
mon = data0['month']
pco2_0 = data0['pco2_sw']
#%%

clim = pco2_0.sel(lat = slice(-30,30),lon = slice(30,120))
#%%
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
    ax.set_title(f"LDEO v2019 pCO₂ EOF mode {i+1} (explains {frac_c[i].values*100:.2f}% variance", pad = 30)

