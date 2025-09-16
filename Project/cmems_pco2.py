#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 13:10:51 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

file = '/home/bobco-08/Desktop/24cl05012/CO2/data/cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m_spco2_40.00E-110.00E_40.00S-30.00N_2022-01-01-2025-09-19.nc'

data = xr.open_dataset(file)

lat =  data['latitude']
lon = data['longitude']
time = data['time']
pco2_c = data['spco2']

pco2_cmems = pco2_c.sel(time = slice('2022-01-01','2023-12-31'))

#%%
months = pco2_cmems['time'].dt.month

# Build season labels as DataArray
seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8]), "JJA", seasons)
seasons = xr.where(months.isin([9, 10, 11]), "SON", seasons)

# Assign as coordinate (must use .data since 'seasons' is a DataArray)
pco2_cmems = pco2_cmems.assign_coords(
    custom_season=("time", seasons.data))

pco2sen_mean = pco2_cmems.groupby('custom_season').mean(dim =('time'))

pco2sen_mean = pco2sen_mean.where(pco2sen_mean != 0)

# define the climatological order you want
season_order = ["DJF", "MAM", "JJA", "SON"]

pco2sen_mean = pco2sen_mean.reindex({"custom_season": season_order})

#%%

fig, axs = plt.subplots(2, 2, figsize=(15,13),dpi = 150 ,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':-0.15,'hspace': 0.19})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    levels = np.linspace(34,45,20)
    im = ax.contourf(pco2_cmems.longitude, pco2_cmems.latitude, pco2sen_mean[i], 
                     levels ,cmap= 'PiYG' ,transform=ccrs.PlateCarree(), extend ='both') #RdYlBu_r nice colormap
    
    contours = ax.contour(pco2_cmems.longitude, pco2_cmems.latitude, pco2sen_mean[i],
                          colors='black', linewidths=0.6, levels= 40)
    ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours
    ax.coastlines()

    ax.add_feature(cf.LAND,linewidth=0.5,color ='grey')
  
 
# Used to add the subplot title inside the plot    
    ax.text(0.75, 0.98, f"Mean ({season[i]})",          # x, y in axes fraction (0–1)
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



plt.suptitle('[cmems]Seasonal climatology of pCO₂ in Indian Ocean (2022-01-01 to 2023-12-31)', fontsize=16, y= 0.95)
plt.show()

#%%

pval  = pco2sen_mean.values

p5, p95 = np.nanpercentile(pval, [5, 95])  # ignores NaNs
print(p5, p95)