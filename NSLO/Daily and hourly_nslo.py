# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:56:34 2025

@author: user
"""


import numpy as np
import xarray as xr

file1 ="C:/Users/user/Documents/24CL05012/nslo/data/data_stream-oper_stepType-instant.nc"
file2 ="C:/Users/user/Documents/24CL05012/nslo/data/u10_daily.nc"
file3 ="C:/Users/user/Documents/24CL05012/nslo/data/v10_daily.nc"


data1 = xr.open_dataset(file1)
data2 = xr.open_dataset(file2)
data3 = xr.open_dataset(file3)



lat = data1['latitude']
lon = data1['longitude']
v = data1['v10']
u = data1['u10']
v1 = data3['v10']
u1 = data2['u10']
time = data1['valid_time']


ws = np.sqrt((u**2)+(v**2))
ws1 = np.sqrt((u1**2)+(v1**2))


latv = lat.values
lonv = lon.values
timev = time.values

import matplotlib.pyplot as plt

ws_aot = ws.sel(latitude = slice(10,-10),longitude = slice(50,70)).mean(dim=('longitude','latitude'))
ws_bobt = ws.sel(latitude = slice(0,-10),longitude = slice(90,110)).mean(dim=('longitude','latitude'))

ws1_aot = ws1.sel(latitude = slice(10,-10),longitude = slice(50,70)).mean(dim=('longitude','latitude'))
ws1_bobt = ws1.sel(latitude = slice(0,-10),longitude = slice(90,110)).mean(dim=('longitude','latitude'))


plt.figure()
plt.plot(ws_aot,label = 'AS', color ='blue')
plt.plot(ws_bobt,label = 'BOB', color ='red')
plt.plot(ws1_aot,label = 'AS', color ='pink')
plt.plot(ws1_bobt,label = 'BOB', color ='green')

plt.xlabel('Time',fontsize = 12)
plt.ylabel('Wind Speed(m/s)',fontsize = 12)
plt.title('Wind speed time series',fontsize = 15)
plt.grid()
plt.legend()

#%%
file1 ="C:/Users/user/Documents/24CL05012/nslo/data/t2m_daily.nc"


data5 = xr.open_dataset(file1)


lat = data5['latitude']
lon = data5['longitude']
temp1 = data5['t2m']-273
time = data5['valid_time']
temp = data1['t2m']-273

t2m_aot = temp.sel(latitude = slice(10,-10),longitude = slice(50,70)).mean(dim=('longitude','latitude'))
t2m_bobt = temp.sel(latitude = slice(0,-10),longitude = slice(90,110)).mean(dim=('longitude','latitude'))

t2m_aot1 = temp1.sel(latitude = slice(10,-10),longitude = slice(50,70)).mean(dim=('longitude','latitude'))
t2m_bobt1 = temp1.sel(latitude = slice(0,-10),longitude = slice(90,110)).mean(dim=('longitude','latitude'))


plt.figure(2, figsize=(12,6),dpi =100)
plt.plot(t2m_aot,label = 'AS', color ='blue')
plt.plot(t2m_bobt,label = 'BOB', color ='red')
plt.plot(t2m_aot1,label = 'AS', color ='blue')
plt.plot(t2m_bobt1,label = 'BOB', color ='red')
plt.xlabel('Time')
plt.ylabel('temperature($^{\circ C}$) ')
plt.title('2-meter temperature time series')
plt.grid()
plt.legend()
