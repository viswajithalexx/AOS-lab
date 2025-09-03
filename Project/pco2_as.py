#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  1 16:29:45 2025

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

pco2_as = pco2.sel(lat = slice(-2,14),lon = slice(40,68))
#%%

pco2_tm = pco2_as.mean(dim = ('lat','lon'))

plt.figure()
plt.plot(pco2_as.mtime,pco2_tm)
plt.grid()
plt.show()
#%%
months = pco2_as['mtime'].dt.month

# Build season labels as DataArray
seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8]), "JJA", seasons)
seasons = xr.where(months.isin([9, 10, 11]), "SON", seasons)

# Assign as coordinate (must use .data since 'seasons' is a DataArray)
pco2_as = pco2_as.assign_coords(
    custom_season=("mtime", seasons.data))

pco2sen_mean = pco2_as.groupby('custom_season').mean()

pco2sen_mean = pco2sen_mean.where(pco2sen_mean != 0)

# define the climatological order you want
season_order = ["DJF", "MAM", "JJA", "SON"]

pco2sen_mean = pco2sen_mean.reindex({"custom_season": season_order})

#%%

fig, axs = plt.subplots(2, 2, figsize=(15,13), subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':0.15,'hspace': 0.19})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    levels = np.linspace(350,400,20)
    im = ax.contourf(pco2_as.lon, pco2_as.lat, pco2sen_mean[i], 
                     levels,cmap= 'rainbow' ,transform=ccrs.PlateCarree(), extend ='both')
    
    contours = ax.contour(pco2_as.lon, pco2_as.lat, pco2sen_mean[i],
                          colors='black', linewidths=0.6, levels= 20)
    ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours
    ax.coastlines()
    ax.add_feature(cf.BORDERS, linewidth=0.5 )
    ax.add_feature(cf.LAND,color ='grey', zorder =15)
 
# Used to add the subplot title inside the plot    
    ax.text(0.75, 0.98, season[i],          # x, y in axes fraction (0–1)
        transform=ax.transAxes,         # interpret coords relative to axes
        fontsize=12, zorder= 18,fontweight='bold',
        va='top', ha='right',)            # vertical & horizontal alignment
    gl = ax.gridlines(draw_labels = True,linewidth = 0.5 , color = 'grey', alpha =0.5)
    gl.right_labels = False
    gl.top_labels = False
    
    # remove *bottom labels* only for first row (i = 0,1)
    if i in [0, 1]:
        gl.bottom_labels = False
    if i in[1,3]:
        gl.left_labels = False
        

    


# Add colorbar
cbar_ax = fig.add_axes([0.25, 0.50, 0.5, 0.015]) # left,bottom,width,height
cbar = fig.colorbar(im, cax=cbar_ax,orientation = 'horizontal', label= '\u00b5atm')



plt.suptitle('Seasonal climatology of pCO₂  in AS (1957–2023)', fontsize=16, y= 0.95)
plt.show()
