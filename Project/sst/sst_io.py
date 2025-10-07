#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 18:15:52 2025

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

file = '/home/bobco-08/Desktop/24cl05012/CO2/data/sst_mon_io_19402025.nc'

sst_ds = xr.open_dataset(file)

lat = sst_ds['latitude']
lon = sst_ds['longitude']
sst = sst_ds['sst']-273.15
time = sst_ds['valid_time']

sst_io = sst.sel(latitude = slice(30,0),longitude = slice(38,110),valid_time = slice('1980-01-01','2019-12-31'))

#%% sst indian ocean

sst_temporal = sst_io.mean(dim = ('latitude','longitude'))
sst_io_mon = sst_temporal.groupby('valid_time.month').mean()
month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

plt.figure(figsize=(10, 8), dpi=150)
plt.plot(month, sst_io_mon, marker='o', linestyle='-', color='r')  
plt.ylim(26,31)
plt.xlabel("Month")
plt.ylabel("SST (°C)")
plt.title("Climatology of monthly SST over Indian Ocean(0°N–30°N, 30°E–110°E)")
plt.grid(True)
plt.show()

#%% sst arabian sea

sst_as = sst.sel(latitude = slice(30,0),longitude = slice(38,78),valid_time = slice('1980-01-01','2019-12-31'))

sst_temporal_as = sst_as.mean(dim = ('latitude','longitude'))
sst_as_mon = sst_temporal_as.groupby('valid_time.month').mean()
month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

plt.figure(figsize=(10, 8), dpi=150)
plt.plot(month, sst_as_mon, marker='o', linestyle='-', color='r') 
plt.ylim(26,31) 
plt.xlabel("Month")
plt.ylabel("SST (°C)")
plt.title("Climatology of monthly SST over Arabian Sea (0°N–30°N, 30°E–78°E)")
plt.grid(True)
plt.show()

#%% sst bay of bengal

sst_bob = sst.sel(latitude = slice(30,0),longitude = slice(78,110),valid_time = slice('1980-01-01','2019-12-31'))

sst_temporal_bob = sst_bob.mean(dim = ('latitude','longitude'))
sst_bob_mon = sst_temporal_bob.groupby('valid_time.month').mean()
month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

plt.figure(figsize=(10, 8), dpi=150)
plt.plot(month, sst_bob_mon, marker='o', linestyle='-', color='r')  
plt.ylim(26,31)
plt.xlabel("Month")
plt.ylabel("SST (°C)")
plt.title("Climatology of monthly SST over Bay of Bengal(0°N–30°N, 78°E–110°E)")
plt.grid(True)
plt.show()
