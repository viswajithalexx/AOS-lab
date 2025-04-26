# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 09:55:16 2025

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr


file = "C:/Users/user/Documents/24CL05012/nsla/olr-daily.nc"

data = xr.open_dataset(file)

lon = data['lon']
lat =  data['lat']
time = data['time']
olr =  data['olr']

lonv =lon.values
latv =lat.values
# timev =time.values


olr_lat = olr.sel(time=olr['time.month']<=3,lat = slice(-15,0)).mean(dim=('lat'))
time_olr = time.sel(time= time['time.month']<=3)
plt.figure(figsize=(12,8),dpi =100)
plt.pcolor(lon,time_olr,olr_lat,cmap = 'jet')
plt.colorbar(label='W/m^2')
plt.xlabel('Longitude')
plt.ylabel('Time')

olr_lon = olr.sel(time= slice('2010-06-01','2010-09-30'),lat =slice(-5,45),lon = slice(75,85)).mean(dim=('lon'))
time_olr2 = time.sel(time= slice('2010-06-01','2010-09-30'))
lat_olr = lat.sel( lat = slice(-5,45))
plt.figure(figsize=(12,8),dpi =100)
plt.pcolor(lat_olr,time_olr2,olr_lon,cmap = 'jet')
plt.colorbar(label='W/m^2')
plt.xlabel('Latitude')
plt.ylabel('Time')

