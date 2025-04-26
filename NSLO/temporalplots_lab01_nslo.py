# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 16:09:14 2025

@author: user
"""

import numpy as np
import xarray as xr

file1 ="C:/Users/user/Documents/24CL05012/nslo/data/u10_daily.nc"
file2 ="C:/Users/user/Documents/24CL05012/nslo/data/v10_daily.nc"
data1 = xr.open_dataset(file1)
data2 = xr.open_dataset(file2)

lat = data1['latitude']
lon = data1['longitude']
v = data2['v10']
u = data1['u10']
time = data1['valid_time']

ws = np.sqrt((u**2)+(v**2))

#wsv = ws.values
latv = lat.values
lonv = lon.values

import matplotlib.pyplot as plt

ws_sio = ws.sel(latitude = slice(20,0),longitude = slice(40,100)).mean(dim=('longitude','latitude'))
# t2m_bobt = temp.sel(latitude = slice(0,-10),longitude = slice(90,110)).mean(dim=('longitude','latitude'))

plt.figure(1)
plt.plot(ws_sio)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Wind speed spatial plot')

#%% time series of t2m

file1 ="C:/Users/user/Documents/24CL05012/nslo/data/t2m_daily.nc"

data3 = xr.open_dataset(file1)


lat = data3['latitude']
lon = data3['longitude']
temp = data3['t2m']-273
time = data3['valid_time']


latv = lat.values
lonv = lon.values

import matplotlib.pyplot as plt

t2m_aot = temp.sel(latitude = slice(10,-10),longitude = slice(50,70)).mean(dim=('longitude','latitude'))
t2m_bobt = temp.sel(latitude = slice(0,-10),longitude = slice(90,110)).mean(dim=('longitude','latitude'))


plt.figure(1, figsize=(12,6),dpi =100)
plt.plot(t2m_aot,label = 'AO', color ='blue')
plt.plot(t2m_bobt,label = 'BOB', color ='red')
plt.xlabel('Time')
plt.ylabel('temperature(')
plt.title('2-meter temperature time series')
plt.grid()
plt.legend()
#%% time series of sst
file1 ="C:/Users/user/Documents/24CL05012/nslo/data/sst-daily.nc"

data4 = xr.open_dataset(file1)


lat = data4['latitude']
lon = data4['longitude']
sst_temp = data4['sst']-273
time = data4['valid_time']


latv = lat.values
lonv = lon.values

import matplotlib.pyplot as plt

sst_aot = sst_temp.sel(latitude = slice(10,-10),longitude = slice(50,70)).mean(dim=('longitude','latitude'))
sst_bobt = sst_temp.sel(latitude = slice(0,-10),longitude = slice(90,110)).mean(dim=('longitude','latitude'))


plt.figure(1, figsize=(12,6),dpi =100)
plt.plot(sst_aot,label = 'AO', color ='blue')
plt.plot(sst_bobt,label = 'BOB', color ='red')
plt.xlabel('Time')
plt.ylabel('temperature')
plt.title('sea surface temperature time series')
plt.grid()
plt.legend()
#%% time series of sp
file1 ="C:/Users/user/Documents/24CL05012/nslo/data/sp_daily.nc"

data5 = xr.open_dataset(file1)


lat = data5['latitude']
lon = data5['longitude']
sp = data5['sp']/100
time = data5['valid_time']


latv = lat.values
lonv = lon.values

import matplotlib.pyplot as plt

sp_aot = sp.sel(latitude = slice(10,-10),longitude = slice(50,70)).mean(dim=('longitude','latitude'))
sp_bobt = sp.sel(latitude = slice(0,-10),longitude = slice(90,110)).mean(dim=('longitude','latitude'))


plt.figure(1, figsize=(12,6),dpi =100)
plt.plot(sp_aot,label = 'AO', color ='blue')
plt.plot(sp_bobt,label = 'BOB', color ='red')
plt.xlabel('Time')
plt.ylabel('Surface Pressure (hPa)')
plt.title('Surface Pressure time series')
plt.grid()
plt.legend()

