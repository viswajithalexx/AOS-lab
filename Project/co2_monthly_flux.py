#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 17:35:35 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature



file = "/home/bobco-08/Desktop/24cl05012/CO2/data/oc_v2024E.flux.nc"
file2 = "/home/bobco-08/Desktop/24cl05012/CO2/data/monthly_sst_io.nc"

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
area = data['dxyp']


# co2_io = co2_ocean.sel(lat = slice(-40,30),lon = slice(30,120),mtime = slice("2000-1-1","2023-12-31")).mean(dim=('mtime'))

# lon_io = lon.sel(lon =slice(30,120))
# lat_io = lat.sel(lat = slice(-40,30))

co2_flux = co2_ocean/area

co2_flux = co2_ocean.sel(lat = slice(-40,30),lon = slice(40,120))*(10**15)

# co2_flux = co2_flux/12.01


co2mon = co2_flux.groupby('mtime.month').mean(dim = 'mtime')
co2mon = co2mon.where(co2mon != 0)


#%%

fig, axs = plt.subplots(3, 4, figsize=(15,14), subplot_kw={'projection': ccrs.PlateCarree()})
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

for i, ax in enumerate(axs.flat):
    
    im = ax.contourf(co2_flux.lon, co2_flux.lat, co2mon[i], 
                     levels=20 , cmap='rainbow', transform=ccrs.PlateCarree())
    
    ax.set_title(months[i])
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.set_xticks([])
    ax.set_yticks([])

# Add colorbar
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5]) # left,bottom,width,height
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label(label =r'g C m$^{-2}$ yr$^{-1}$',labelpad = 10)


plt.suptitle('Monthly Mean CO₂ Flux (1957–2023)', fontsize=16)
plt.tight_layout(rect=[0, 0, 0.9, 0.95])
plt.show()