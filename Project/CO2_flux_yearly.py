# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 18:11:30 2025

@author: user
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature

file = "C:/Users/user/Documents/24CL05012/2nd Year/Project/code/data/oc_v2024E.flux_yearly_3Bands30.nc"

data = xr.open_dataset(file)

# variables 
co2flux = data['co2flux']
lon = data['lon']
lat = data['lat']
reg = data['reg']
time = data['mtime']

co2v = co2flux.values

flux1 = co2flux.sel(reg = 4) #30N - 90N
flux2 = co2flux.sel(reg = 5) #30S -30N
flux3 = co2flux.sel(reg = 6) #90S - 30S

#%%

# yearly co2 flux 1957-07-01 to 2023-07-01
region_names = {
    4: "30°N - 90°N",
    5: "30°S -30°N",
    6: "90°S - 30°S"

}

for i in range(4,7):
    
     co2_region = co2flux.sel(reg = i)
     plt.figure()
     plt.plot(time,co2_region, label = region_names[i])
     plt.xlabel('Years')
     plt.ylabel('PgC/yr')
     plt.title(f"1957 -2023 yearly CO2 flux: {region_names[i]}")
     plt.grid(True)
#%%

         
co2_region = co2flux.sel(reg = 8)
plt.figure()
plt.plot(time,co2_region, label = 'Global')
plt.xlabel('Years')
plt.ylabel('PgC/yr')
plt.title("1957 -2023 yearly CO2 flux: Global")
plt.grid(True)
  
#%%

plt.figure()
plt.plot(time,flux1, label ='30N-90N')
plt.plot(time,flux2,label ='30S-30N')
plt.plot(time,flux3,label ='90S-30S')
plt.xlabel('Years')
plt.ylabel('PgC/yr')
plt.grid()
plt.legend()
#%%

# without tropics

plt.figure()
plt.plot(time,flux1, label ='30N-90N')
plt.plot(time,flux3,label ='90S-30S')
plt.xlabel('Years')
plt.ylabel('PgC/yr')
plt.title('without tropics')
plt.grid()
plt.legend()

