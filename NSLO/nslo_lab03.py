# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 18:19:03 2025

@author: user
"""
import numpy as np


def chlorophyll_value(chl,z):
    chl_midpoint =[]
    depth_z =[]
    chlorophyll_sum =[]

    for i in range(len(chl)-1):
             chl_m = ((chl[i+1]+chl[i])/2)
             chl_midpoint.append(chl_m)
             depth = z[i+1]-z[i]
             depth_z.append(depth)
 
    for i in range(len(chl_midpoint)):
           chlorophyll = (chl_midpoint[i]*depth_z[i])
           chlorophyll_sum.append(chlorophyll) 
    chlorophyll_value = np.sum(chlorophyll_sum)
    return chlorophyll_value


chl = np.array([0.35,0.40,0.50,0.55,0.48,0.30]) 
z = np.array([5,10,20,40,70,100])

chl_v = chlorophyll_value(chl,z)
print('value of chlorophyll Concentration=',chl_v)

#%%
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

file = "C:/Users/user/Documents/24CL05012/nslo/data/d.nc"
data = xr.open_dataset(file)

lon = data['longitude']
lat = data['latitude']
time = data['valid_time']
sst  = data['sst']
sst_int = []
latv =lat.values
lat_sub = np.arange(8,14,0.25)
for i in np.arange(8,14,0.25):
    
      
        sst_selection =sst.sel(latitude = i,longitude = slice(82, 90)).mean(dim=('valid_time'))
        sst_values = sst_selection.values-273.15
        longitude_values = np.linspace(82,90,len(sst_values))
        chl = sst_values
        z = longitude_values
        chl_v = chlorophyll_value(chl,z)
        sst_int.append(chl_v)
print('sst integrated along latitudes',np.mean(sst_int),'°C')       
plt.figure(figsize =(12,8),dpi=100)
plt.plot(lat_sub,sst_int)
plt.grid()
plt.xlabel('Latitude (°E)', fontsize=12)
plt.ylabel('Sea Surface Temperature integrated °C ')
plt.title('Sea surface temperature Integrated along Latitudes')

#%%
sst_int_lon = []
latv =lat.values
lon_sub = np.arange(82,90,0.25)

for i in np.arange(82,90,0.25):

        sst_selection =sst.sel(latitude = slice(13,8),longitude = i).mean(dim=('valid_time'))
        sst_values = sst_selection.values-273.15
        latitude_values = np.linspace(8,13,len(sst_values))
        chl = sst_values
        z = latitude_values
        chl_v = chlorophyll_value(chl,z)
        sst_int_lon.append(chl_v)
print('sst integrated along longitudes',np.mean(sst_int_lon),'°C')
plt.figure(figsize =(12,8),dpi=100)
plt.plot(lon_sub,sst_int_lon)
plt.grid()
plt.xlabel('Longitude (°E)', fontsize=12)
plt.ylabel('Sea Surface Temperature °C')
plt.title('Sea surface temperature integrated along Longitude')
#%%
