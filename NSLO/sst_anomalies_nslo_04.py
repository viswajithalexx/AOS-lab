# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 11:36:49 2025

@author: HP
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

path = "C:/Users/user/Documents/24CL05012/nslo/data/hadsst_global_1990_2023.nc"
data = xr.open_dataset(path)

lon = data['LON']
lat = data['LAT']
time = data['TIME']
sst = data['SST']
lonv = lon.values
latv = lat.values
timev= time.values
sst_pacific = sst.sel(LON = slice(-170,-120),LAT = slice(-5,5))

sst_clim = sst_pacific.groupby('TIME.month').mean()
sst_monthly = sst_pacific.groupby('TIME.month')

sst_anomaly = sst_monthly-sst_clim

sst_wdetrend = sst_anomaly.mean(dim=('LAT','LON'))

plt.figure(figsize=(18,10), dpi=100)
x= timev
y= sst_wdetrend
plt.plot(x,y,color='green')
plt.axhline(0.5,linestyle ='--',color ='r')
plt.axhline(0,linestyle ='-',color ='black')
plt.axhline(-0.5,linestyle ='--',color ='b')
plt.fill_between(x, y, 0.5, where=(y > 0.5), color='red', alpha=0.5)
plt.fill_between(x, y, -0.5, where=(y <-0.5), color='b', alpha=0.5)
plt.xlabel('YEARS')
plt.ylabel('DETREND SST')
plt.title('SST ANOMALIES WITHOUT DETREND')

#%% sst  detrend 

cof = sst_pacific.polyfit(dim='TIME', deg=1)
trends = xr.polyval(sst_pacific['TIME'], cof.polyfit_coefficients)
dtrend=sst_pacific-trends

sst_clim=dtrend.groupby('TIME.month').mean()

sst_monthly = dtrend.groupby('TIME.month')

sst_anomalies =sst_monthly-sst_clim

sst_detrend = sst_anomalies.mean(dim =('LAT','LON'))

plt.figure(figsize=(18,10), dpi=100)
x= timev
y= sst_detrend
plt.plot(x,y,color='black')
plt.axhline(0.5,linestyle ='--',color ='r')
plt.axhline(0,linestyle ='-',color ='black')
plt.axhline(-0.5,linestyle ='--',color ='b')
plt.fill_between(x, y, 0.5, where=(y > 0.5), color='red', alpha=0.5)
plt.fill_between(x, y, -0.5, where=(y <-0.5), color='b', alpha=0.5)
plt.xlabel('YEARS')
plt.ylabel('DETREND SST')
plt.title('SST ANOMALIES WITH DETREND')
#%%
plt.figure(figsize=(18,10), dpi=100)
plt.plot(timev,sst_detrend,color='red',label ='with detrend')
plt.plot(timev,sst_wdetrend,color='blue',label= 'without detrend')
plt.axhline(0.5,linestyle ='--',color ='r')
plt.axhline(0,linestyle ='-',color ='black')
plt.axhline(-0.5,linestyle ='--',color ='b')
plt.xlabel('YEARS')
plt.ylabel('DETREND SST')
plt.title('SST ANOMALIES WITH DETREND vs WITHOUT DETREND')
plt.legend()