# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 11:06:25 2025

@author: user
"""

import numpy as np
import xarray as xr
file ="C:/Users/user/Documents/24CL05012/nsla/WIND_DATA/U_V_Wind_June_1979_2021_ERA5_Hr_BOB.nc"

data = xr.open_dataset(file)

lon = data['longitude']
lat = data['latitude']
time = data['time']
u = data['u10']
v = data['v10']

lonv =lon.values
latv =lat.values
timev = time.values

u10 = u.sel(time = slice('1979-06-01','1979-06-10'))
v10 = v.sel(time = slice('1979-06-01','1979-06-10'))

dx =  (lonv[1]-lonv[0])* 110 * 1000
dy =  (latv[1]-latv[0])* 110 * 1000 

#%%

du_dx = np.gradient(u10,axis=2)/dx
dv_dy = np.gradient(v10,axis=1)/dy
du_dx_mean = np.mean(du_dx,axis =0)
dv_dy_mean = np.mean(dv_dy,axis=0)

div = du_dx_mean+dv_dy_mean
du_dy = np.gradient(u10,axis=2)/dy
dv_dx = np.gradient(v10,axis=1)/dx
vor = dv_dx-du_dy
vor = np.mean(vor,axis=0)

#%%

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
ops =[du_dx_mean,dv_dy_mean,div,vor]
fig,wind_subplots = plt.subplots(nrows=2,ncols=2,figsize=(12,8),subplot_kw={'projection':ccrs.PlateCarree()})
wind_subplots = wind_subplots.flatten()
for i in range(0,3):
    ax =wind_subplots[i]  
    wind = contourf(du_dx_mean,cmap='jet',add_colorbar=False, levels=22, extend="both",zorder=5)
    ax.add_feature(cf.LAND,color="black",zorder=11)    
plt.colorbar(wind, ax=wind_subplots, orientation="vertical", pad=0.05)   
#%%
 