# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 14:10:56 2025

@author: HP
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cf


file= "C:/Users/user/Documents/24CL05012/nslo/data/sst-daily.nc"

data = xr.open_dataset(file)

lon = data['longitude']
lat = data['latitude']
time = data['valid_time']
sst  = data['sst']

# day1 = time[0]

latv =lat.values
sstv = sst.values

fig,sst_subplots=plt.subplots(5,2,figsize=(20,18),subplot_kw={'projection':ccrs.Mercator()})
sst_subplots = sst_subplots.flatten()
plt.subplots_adjust(hspace= 0.9)

for i in range(0,10): 
    sst_spatial = sst.sel(valid_time=time[i], latitude=slice(30, 5), longitude=slice(50, 110))
    sstd_s = sst_spatial-273.15
    ax = sst_subplots[i]
    contour = sstd_s.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap='jet',vmax =36,vmin =18,add_colorbar=False, levels=22, extend="both",zorder=5) #extend  makes values below vmin one colour and above one colour
    ax.add_feature(cf.LAND,color="black",zorder=11)
    ax.set_title(f"Day{i+1}")
    gl =ax.gridlines(draw_labels =True)
plt.colorbar(contour, ax=sst_subplots, orientation="vertical", pad=0.05, label="Temperature (°C)")
plt.savefig('sst_31.png')

 #%%
fig1,sst_subplots1=plt.subplots(5,2,figsize=(20,18),subplot_kw={'projection':ccrs.Mercator()})
sst_subplots1 = sst_subplots1.flatten()
plt.subplots_adjust(hspace = 0.9)

for i in range(10,20): 
     sst_spatial1 = sst.sel(valid_time=time[i], latitude=slice(30, 5), longitude=slice(50, 120))
     sstd_s1 = sst_spatial1-273.15
     ax1 = sst_subplots1[i-10]
     ax1.gridlines(draw_labels =True)
     contour1 = sstd_s1.plot.contourf(ax=ax1, transform=ccrs.PlateCarree(), cmap='jet',vmax =36,vmin =18,add_colorbar=False, levels=22, extend="both",zorder=5)
     ax1.set_title(f"Day{i+1}")
     ax.add_feature(cf.LAND,color="black",zorder=11)
     ax1.add_feature(cf.LAND,color="black",zorder=11)
plt.colorbar(contour1, ax=sst_subplots1, orientation="vertical", pad=0.05, label="Temperature (°C)")
plt.savefig('sst_32.png')
#%%
fig2,sst_subplots2=plt.subplots(5,2,figsize=(20,18),subplot_kw={'projection':ccrs.Mercator()})
sst_subplots2 = sst_subplots2.flatten()
plt.subplots_adjust(hspace = 0.9)

for i in range(20,30): 
     sst_spatial2 = sst.sel(valid_time=time[i], latitude=slice(30, 5), longitude=slice(50, 120))
     sstd_s2 = sst_spatial2-273.15
     ax2 = sst_subplots2[i-20]
     ax2.gridlines(draw_labels =True)
     contour2 = sstd_s2.plot.contourf(ax=ax2, transform=ccrs.PlateCarree(), cmap='jet',vmax =36,vmin =18,add_colorbar=False, levels=22, extend="both",zorder=5)
     ax2.add_feature(cf.LAND,color="black",zorder=11)
     ax2.set_title(f"Day{i+1}")

plt.colorbar(contour2, ax=sst_subplots2, orientation="vertical", pad=0.05, label="Temperature (°C)")

plt.savefig('sst_33.png')

 
#%%
   
sst_max =[]
sst_min =[]
sstd_min_july = []
sstd_max_july = []

for n in range(0,31):
    sst_time = sst.sel(valid_time = time[n],latitude = slice(13, 8),longitude = slice(82, 90))
    sstr_c = sst_time-273.15
    sstd_min_july = np.min(sstr_c).values
    sstd_max_july = np.max(sstr_c).values
    sst_max.append(sstd_max_july)
    sst_min.append(sstd_min_july)
    
    
  
plt.figure(figsize=(15,10),dpi=100)
plt.plot(time,sst_max, label = 'Daily Highest SST ', color = 'blue',marker = 'o')  
plt.plot(time,sst_min, label = 'Daily Lowest SST ', color = 'green',marker = 'x')
plt.xlabel('Dates of july')
plt.ylabel('SST (in degree Celsius)')
plt.title('Sea Surface Temperature both Daily Maximum and Minimum for July, 2023')
plt.savefig('sst_max_min.png')
plt.legend()
plt.show() 