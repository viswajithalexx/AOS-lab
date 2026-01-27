#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 12:55:55 2026

@author: bobco-08
"""

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt

data1 = xr.open_dataset('/home/bobco-08/24cl05012/TDcv/rhum.2001.nc')
data2 = xr.open_dataset('/home/bobco-08/24cl05012/TDcv/rhum.2002.nc')

# --------------------
# Select 850 hPa level
# --------------------
rh1_850 = data1['rhum'].sel(level=850)
rh2_850 = data2['rhum'].sel(level=850)

# --------------------
# Select DJF months
# Dec 2001 + Jan Feb 2002
# --------------------
dec_2001 = rh1_850.sel(time=rh1_850['time'].dt.month == 12)
janfeb_2002 = rh2_850.sel(time=rh2_850['time'].dt.month.isin([1, 2]))

# --------------------
# Combine and take mean
# --------------------
djf = xr.concat([dec_2001, janfeb_2002], dim='time')
djf_mean = djf.mean(dim='time')

# --------------------
# Plot spatial map
# --------------------
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

pcm = ax.pcolormesh(
    djf_mean['lon'],
    djf_mean['lat'],
    djf_mean,
    shading='auto',
    cmap='YlGnBu',
    transform=ccrs.PlateCarree()
)

plt.colorbar(pcm, ax=ax, label='Relative Humidity (%)')
ax.coastlines()
ax.add_feature(cf.BORDERS, linewidth=0.5)
ax.set_title('DJF Mean Relative Humidity at 850 hPa (Dec 2001 – Feb 2002)')

plt.show()


#%%

rh1 = data1['rhum']
rh2 = data2['rhum']

dec_2001 = rh1.sel(time=rh1['time'].dt.month == 12)
janfeb_2002 = rh2.sel(time=rh2['time'].dt.month.isin([1, 2]))

djf_pl = xr.concat([dec_2001, janfeb_2002], dim='time')

# --------------------
# Mean over time (DJF mean)
# --------------------
djf_mean = djf_pl.mean(dim='time')

# --------------------
# Zonal mean (average over longitude)
# --------------------
lat_pres = djf_mean.mean(dim='lon')

# --------------------
# Plot latitude–pressure cross section
# --------------------
plt.figure(figsize=(10, 6))

# Filled contours
cf = plt.contourf(
    lat_pres['lat'],
    lat_pres['level'],
    lat_pres,
    levels=20,
    cmap='YlGnBu',
    extend='both'
)

# Contour lines
cs = plt.contour(
    lat_pres['lat'],
    lat_pres['level'],
    lat_pres,
    levels=10,
    colors='k',
    linewidths=0.7
)

plt.clabel(cs, fmt='%d', fontsize=9)   # label contour lines

plt.gca().invert_yaxis()
plt.colorbar(cf, label='Relative Humidity (%)')

plt.xlabel('Latitude')
plt.ylabel('Pressure (hPa)')
plt.title('Latitude–Pressure Cross Section of RH (DJF: Dec 2001 – Feb 2002)')

plt.show()

