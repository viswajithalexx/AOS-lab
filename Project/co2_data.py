# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 16:21:18 2025

@author: HP
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature


file = "/home/bobco-08/Desktop/24cl05012/CO2/codes/oc_v2024E.flux.nc"
file2 = "/home/bobco-08/Desktop/24cl05012/CO2/codes/monthly_sst_io.nc"

data = xr.open_dataset(file,decode_timedelta=True)
data2 = xr.open_dataset(file2,decode_timedelta = True)

#SST data variables

lon_sst = data2['longitude']
lat_sst = data2['latitude']
time_sst = data2['valid_time']
sst  = data2['sst']


#CO2 data variables

lat =  data['lat']
lon = data['lon']
time = data['mtime']
co2_ocean = data['co2flux_ocean']

co2_io = co2_ocean.sel(lat = slice(-40,30),lon = slice(30,120),mtime = slice("2000-1-1","2023-12-31")).mean(dim=('mtime'))

lon_io = lon.sel(lon =slice(30,120))
lat_io = lat.sel(lat = slice(-40,30))

#%% spatial

# Create figure and axis with PlateCarree projection
x = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot coastlines and land
ax.coastlines()
ax.add_feature(cfeature.LAND, color='lightgray', zorder=11)

# Plot the data with 'coolwarm' colormap
im = ax.contourf(lon_io, lat_io, co2_io, transform=ccrs.PlateCarree())

# Add gridlines with labels
gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10}
gl.ylabel_style = {'size': 10}

# Add colorbar with label
cbar = plt.colorbar(im, orientation='vertical', pad=0.05, aspect=30)
cbar.set_label("CO₂ Flux (PgC/yr)", fontsize=12)

# Show plot
plt.title("CO₂ Flux over Indian Ocean", fontsize=14)
plt.tight_layout()
plt.show()

#%% IOD and CO2 flux

co2_flux = co2_ocean.sel(lat = slice(-40,30),lon = slice(30,120))

# Remove all Feb 29s from the time axis
co2_flux = co2_flux.sel(mtime=~((co2_flux['mtime.month'] == 2) & (co2_flux['mtime.day'] == 29)))


co2_monthly = co2_flux.resample(mtime = '1MS').mean()

co2_clim = co2_monthly.groupby('mtime.month').mean('mtime')

co2_anomaly = co2_monthly.groupby('mtime.month') - co2_clim

co2a_t = co2_anomaly.mean(dim=('lat','lon'))


#IOD index

sst_modified = sst.sel(valid_time = slice('1957-1-1','2023-12-31'))
sst_clim = sst_modified.groupby('valid_time.month').mean()

sst_anomaly = sst_modified.groupby('valid_time.month')- sst_clim

wbox = sst_anomaly.sel(latitude = slice(10,-10),longitude = slice(50,70)).mean(dim = ('latitude','longitude'))
ebox = sst_anomaly.sel(latitude =slice(0,-10),longitude = slice(90,110)).mean(dim = ('latitude','longitude'))

iod = wbox - ebox


#%% Plotting IOD and CO2
iod_time = np.linspace(1957,2023,804)

plt.figure()
plt.plot(iod_time,iod)
plt.plot(iod_time,co2a_t)
plt.show()


