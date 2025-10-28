#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 21:21:42 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
from eofs.xarray import Eof
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean.cm as cm

file2 = '/home/bobco-08/24cl05012/CO2/data/pCO2-Corrected_INCOIS-BIO-ROMS_v2.nc'
data2 = xr.open_dataset(file2)


lat2 = data2['LAT']
lon2 = data2['LON']
time2 = data2['TIME']
pco2_2 = data2['pCO2_Original']


incois = pco2_2.sel(LAT = slice(0,30),LON = slice(30,120))
#%%
di = incois
#%%


di_rm = di.rolling(TIME=3, center=True).mean()
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
# di_detrended.name = 'pCO2_detrended'
# di_detrended.to_netcdf('/home/bobco-08/24cl05012/CO2/data/pco2_incois_detrended.nc')
#%%

detrend = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/pco2_incois_detrended.nc')

di_detrended = detrend['pCO2_detrended']



#%%
di_detrended = di_detrended.interpolate_na(dim= 'TIME',method = 'linear', fill_value="extrapolate")

di_ocean = di_detrended.where(~di_detrended.isnull().all(dim = ('TIME')))



#%%
solver_i = Eof(di_ocean)

# --- Extract the first 3 EOFs and PCs ---
eofs_i  = solver_i.eofsAsCorrelation(neofs=3)   # spatial patterns
pcs_i   = solver_i.pcs(npcs=3, pcscaling=1)     # time series
frac_i  = solver_i.varianceFraction(neigs=3)    # variance explained

# --- Print explained variance ---
for i in range(3):
    print(f"INCOIS pCO₂ EOF mode {i+1} (explains {frac_i[i].values*100:.2f}% variance)"
)
#%% eof spatial

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
