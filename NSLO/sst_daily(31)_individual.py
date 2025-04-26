# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 18:48:44 2025

@author: user
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

day1 = time[0]

latv =lat.values
sstv = sst.values



for i in range(0,31):
    
    sst_spatial = sst.sel(valid_time=time[i], latitude=slice(30, 5), longitude=slice(50, 100))
    sstd_s = sst_spatial-273.15
    plt.figure(figsize=(12,6),dpi=100)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    contour = sstd_s.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap='RdBu_r',vmax =36,vmin =18,add_colorbar=False, levels=22, extend="both",zorder=5)
    ax.add_feature(cf.LAND,color="black",zorder=11)
    ax.set_title(f'Day{i+1}')
    ax.gridlines(draw_labels = False)
    plt.colorbar(contour,orientation="vertical", pad=0.05, label="Temperature (°C)")
    plt.savefig(f'sst_{i+1}.png')
