# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 18:50:45 2025

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

latv = lat.values
lonv = lon.values

ws = np.sqrt((u**2)+(v**2))

import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature
# Select the region over Bay of Bengal and compute the mean
ws_bob = ws.sel(latitude=slice(20, 5), longitude=slice(83,90)).mean(dim='valid_time')

# Extract the latitudes and longitudes from the selected region
lat_bob = ws_bob.latitude
lon_bob = ws_bob.longitude

# Get the wind components for quiver plot
uv = u.sel(latitude=slice(20, 5), longitude=slice(75, 100)).values
vv = v.sel(latitude=slice(20, 5), longitude=slice(75, 100)).values

# Create the meshgrid for contour plot
lon_grid, lat_grid = np.meshgrid(lon_bob, lat_bob)

# Set up the xtick and ytick values based on the latitudes and longitudes
xtick = lon_bob.values
ytick = lat_bob.values

# Plot
plt.figure(figsize=(15, 10), dpi=100)
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.COASTLINE)
ax1.add_feature(cfeature.BORDERS)
ax1.add_feature(cfeature.LAND)

# Use contourf to plot the wind speed
im = ax1.contourf(lon_grid, lat_grid, ws_bob, cmap='jet')

# Use quiver to plot the wind vectors
sp = 3  # sampling for quiver
ax1.quiver(xtick[::sp], ytick[::sp], uv[0, ::sp, ::sp], vv[0, ::sp, ::sp], scale=200, width=0.002)

# Add gridlines and colorbar
ax1.gridlines(visible=True, draw_labels=True)
plt.colorbar(im, label='Wind Speed (m/s)')

# Add a title
ax1.set_title('Wind Speed and Vector Field (Bay of Bengal)')

# Show the plot
plt.show()
