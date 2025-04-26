# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 01:00:46 2025

@author: HP
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

file1 = 'E:/Work/IITBBS/Mr.robot/Python/nslo/data/data_stream-oper_stepType-instant.nc'

data1 = xr.open_dataset(file1)

lat = data1['latitude']
lon = data1['longitude']
v = data1['v10']
u = data1['u10']
sp = data1['sp'] / 100
sst = data1['sst'] - 273
temp = data1['t2m'] - 273
time = data1['valid_time']

ws = np.sqrt((u ** 2) + (v ** 2))

# Calculate mean values
ws_siom = ws.mean(dim='valid_time')
sst_siom = sst.mean(dim='valid_time')
sp_siom = sp.mean(dim='valid_time')
t2m_siom = temp.mean(dim='valid_time')

fig, axes = plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(15, 10), dpi=100)

# Wind speed
ax = axes[0, 0]
ax.coastlines()
im = ax.contourf(lon, lat, ws_siom)
ax.gridlines(draw_labels=True)
ax.add_feature(cfeature.LAND, color="gray", zorder=11)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Wind Speed (m/s)')
ax.set_title('Mean Wind Speed over Selected Region')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Sea Surface Temperature
ax = axes[0, 1]
ax.coastlines()
im = ax.contourf(lon, lat, sst_siom, cmap='jet')
ax.gridlines(draw_labels=True)
ax.add_feature(cfeature.LAND, color="gray", zorder=11)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Sea Surface Temperature (°C)')
ax.set_title('Sea Surface Temperature over NIO')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Sea Surface Pressure
ax = axes[1, 0]
ax.coastlines()
im = ax.contourf(lon, lat, sp_siom, cmap='gist_earth')
ax.gridlines(draw_labels=True)
ax.add_feature(cfeature.LAND, color="gray", zorder=11)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Sea Surface Pressure (hPa)')
ax.set_title('Surface Pressure over NIO')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# 2-meter Temperature
ax = axes[1, 1]
ax.coastlines()
im = ax.contourf(lon, lat, t2m_siom, cmap='jet')
ax.gridlines(draw_labels=True)
ax.add_feature(cfeature.LAND, color="gray", zorder=11)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('2-meter Temperature (°C)')
ax.set_title('2-meter Temperature over NIO')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

plt.tight_layout()
plt.show()
