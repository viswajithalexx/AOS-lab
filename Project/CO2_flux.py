# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature


file = "C:/Users/HP/Desktop/Project/code/oc_v2024E.flux.nc"

data = xr.open_dataset(file)

lat =  data['lat']
lon = data['lon']
time = data['mtime']
co2_ocean = data['co2flux_ocean']

#%%

co2_ov = co2_ocean.values

co2_o1 = co2_ocean.sel(mtime =slice('2013-01-01','2023-12-31'),lat = slice(-10,30),lon = slice(50,100))
time1 = time.sel(mtime= slice('2013-01-01','2023-12-31'))
lon_io = lon.sel(lon = slice(50,100))
lat_io = lat.sel(lat = slice(-10,30))

co2_

#%% line plot

plt.figure()

co2_om = co2_o1.mean(dim=['lat','lon'])

plt.plot(time1,co2_om)


#%% CO2 flux 2013 - 2023

co2_io = co2_o1.mean(dim =['mtime'])


co2_iov = co2_io.values

plt.figure(figsize=(10, 6))

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.LAND, color="gray", zorder=11)
ax.gridlines(draw_labels=True)

# Choose your preferred colormap here
im = ax.contourf(lon_io, lat_io, co2_io, levels=20, cmap="coolwarm")

# Add colorbar with label
cb = plt.colorbar(im, ax=ax, orientation="vertical", shrink=0.8, pad=0.05)
cb.set_label("CO₂ Flux (PgC/yr)")

plt.title("Ocean CO₂ Flux over IO 2013-2023")
plt.tight_layout()
plt.show()

#%% 
