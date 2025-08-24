#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:50:53 2025

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

file = "/home/bobco-08/Desktop/24cl05012/CO2/data/oc_v2024E.pCO2.nc"

data = xr.open_dataset(file)

lat =  data['lat']
lon = data['lon']
time = data['mtime']
pco2 = data['pCO2']
area = data['dxyp']

pco2_flux = pco2.sel(lat = slice(-40,30),lon = slice(40,120))

#%%
months = pco2_flux['mtime'].dt.month

# Build season labels as DataArray
seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8]), "JJA", seasons)
seasons = xr.where(months.isin([9, 10, 11]), "SON", seasons)

# Assign as coordinate (must use .data since 'seasons' is a DataArray)
pco2_flux = pco2_flux.assign_coords(
    custom_season=("mtime", seasons.data))

pco2sen_mean = pco2_flux.groupby('custom_season').mean()

pco2sen_mean = pco2sen_mean.where(pco2sen_mean != 0)

# define the climatological order you want
season_order = ["DJF", "MAM", "JJA", "SON"]

pco2sen_mean = pco2sen_mean.reindex({"custom_season": season_order})

#%%

fig, axs = plt.subplots(2, 2, figsize=(15,13), subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':0.25,'hspace': 0.5})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    im = ax.contourf(pco2_flux.lon, pco2_flux.lat, pco2sen_mean[i], 
                     levels =20,cmap= 'rainbow' ,transform=ccrs.PlateCarree(), extend ='both')
    ax.set_title(season[i])
    ax.coastlines()
    ax.add_feature(cf.BORDERS, linewidth=0.5)
    ax.add_feature(cf.LAND,color ='grey',zorder=11)
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)
    gl.top_labels =  False
    gl.right_labels =  False

# Add colorbar
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.25, 0.50, 0.5, 0.015]) # left,bottom,width,height
cbar = fig.colorbar(im, cax=cbar_ax,orientation = 'horizontal', label= r'uatm')



plt.suptitle('Seasonal climatology of CO₂ Flux (1957–2023)', fontsize=16, y= 0.95)
plt.show()

#%% monthly pCO2

pco2mon = pco2_flux.groupby('mtime.month').mean(dim = 'mtime')
pco2mon = pco2mon.where(pco2mon != 0)



fig, axs = plt.subplots(3, 4, figsize=(15,14), subplot_kw={'projection': ccrs.PlateCarree()})
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

for i, ax in enumerate(axs.flat):
    
    im = ax.contourf(pco2_flux.lon, pco2_flux.lat, pco2mon[i], 
                     levels=20 , cmap='rainbow', transform=ccrs.PlateCarree())
    
    ax.set_title(months[i])
    ax.coastlines()
    ax.add_feature(cf.BORDERS, linewidth=0.5)
    ax.add_feature(cf.LAND,linewidth=0.5, zorder =11)
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.5)
    gl.top_labels =  False
    gl.right_labels =  False
    
# Add colorbar
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5]) # left,bottom,width,height
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label(label =r'uatm',labelpad = 10)


plt.suptitle('Monthly Mean CO₂ Flux (1957–2023)', fontsize=16)
plt.show()





