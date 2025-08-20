#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 16:48:29 2025

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


months = co2_flux["mtime"].dt.month

# Build season labels as DataArray
seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8]), "JJA", seasons)
seasons = xr.where(months.isin([9, 10, 11]), "SON", seasons)

# Assign as coordinate (must use .data since 'seasons' is a DataArray)
co2_flux = co2_flux.assign_coords(
    custom_season=("mtime", seasons.data)
)

co2sen_mean = co2_flux.groupby('custom_season').mean()

# define the climatological order you want
season_order = ["DJF", "MAM", "JJA", "SON"]

co2sen_mean = co2sen_mean.reindex({"custom_season": season_order})

co2sen_mean = co2sen_mean.where(co2sen_mean != 0)

#%%



fig, axs = plt.subplots(2, 2, figsize=(15,13), subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':0.25,'hspace': 0.5})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    im = ax.contourf(co2_flux.lon, co2_flux.lat, co2sen_mean[i], 
                     levels =20,cmap= 'rainbow' ,transform=ccrs.PlateCarree(), extend ='both')
    ax.set_title(season[i])
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND,color ='grey',zorder=11)
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)
    gl.top_labels =  False
    gl.right_labels =  False

# Add colorbar
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.25, 0.50, 0.5, 0.015]) # left,bottom,width,height
cbar = fig.colorbar(im, cax=cbar_ax,orientation = 'horizontal', label= r'g C m$^{-2}$ yr$^{-1}$')



plt.suptitle('Seasonal climatology of CO₂ Flux (1957–2023)', fontsize=16, y= 0.95)
plt.show()



#%%

#DJF
ax=plt.figure(1, figsize=(10,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im=ax.contourf(co2_flux.lon, co2_flux.lat, co2sen_mean[0],levels =20,cmap= 'rainbow' ,transform=ccrs.PlateCarree(), extend ='both')          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im,label= r'g C m$^{-2}$ yr$^{-1}$',pad = 0.13)

plt.title('DJF climatology of CO₂ Flux (1957–2023)', pad = 20,fontsize = 16)
plt.grid(True)


#MAM
ax=plt.figure(2,figsize=(10,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im=ax.contourf(co2_flux.lon, co2_flux.lat, co2sen_mean[1],levels =20,cmap= 'rainbow' ,transform=ccrs.PlateCarree(), extend ='both')          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im,label= r'g C m$^{-2}$ yr$^{-1}$',pad =0.13)

plt.title('MAM climatology of CO₂ Flux (1957–2023)',pad = 20, fontsize = 16)
plt.grid(True)


#JJA
ax=plt.figure(3,figsize=(10,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im=ax.contourf(co2_flux.lon, co2_flux.lat, co2sen_mean[2],levels =20,cmap= 'rainbow' ,transform=ccrs.PlateCarree(), extend ='both')          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im,label= r'g C m$^{-2}$ yr$^{-1}$',pad = 0.13)

plt.title('JJA climatology of CO₂ Flux (1957–2023)',pad = 20,fontsize = 16)
plt.grid(True)
#SON
ax=plt.figure(4,figsize=(10,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im=ax.contourf(co2_flux.lon, co2_flux.lat, co2sen_mean[3],levels =20,cmap= 'rainbow' ,transform=ccrs.PlateCarree(), extend ='both')          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im,label= r'g C m$^{-2}$ yr$^{-1}$',pad = 0.13)


plt.title('SON climatology of CO₂ Flux (1957–2023)',pad =20,fontsize = 16)
plt.grid(True)
plt.show()

#%%
