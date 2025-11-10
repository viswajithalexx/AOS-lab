#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 11:45:19 2025

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean.cm as cm
import cmaps

file = '/home/bobco-08/24cl05012/CO2/data/pCO2-Corrected_INCOIS-BIO-ROMS_v2.nc'

data = xr.open_dataset(file)

lat =  data['LAT']
lon = data['LON']
time = data['TIME']
pco2_ibr_ml = data['pCO2_Original']

#%%

# pco2_ibr = pco2_ibr.sel(LAT = slice(-30,30),LON = slice(32.5,110))

months = pco2_ibr_ml['TIME'].dt.month

# Build season labels as DataArray
seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8,9]), "JJAS", seasons)
seasons = xr.where(months.isin([ 10, 11]), "ON", seasons)

# Assign as coordinate (must use .data since 'seasons' is a DataArray)
pco2_ml = pco2_ibr_ml.assign_coords(
    custom_season=("TIME", seasons.data))

pco2sen_ml = pco2_ml.groupby('custom_season')

pco2sen_ml = pco2sen_ml.mean(dim = ('TIME'))


# define the climatological order you want
season_order = ["DJF", "MAM", "JJAS", "ON"]

pco2sen_ml = pco2sen_ml.reindex({"custom_season": season_order})

#%%

pco2sen_ml.to_netcdf('/home/bobco-08/24cl05012/CO2/data/pCO2-Corrected_INCOIS-BIO-ROMS_seasonal.nc')

#%%

fig, axs = plt.subplots(2, 2, figsize=(15,10),dpi = 300 ,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':-0.002,'hspace': 0.06})
season = ['DJF','MAM','JJAS',"SON"]

for i, ax in enumerate(axs.flat):
    levels = np.arange(300,460,10)
    im = ax.contourf(pco2_ml.LON, pco2_ml.LAT, pco2sen_ml[i], 
                     levels = levels,cmap= cmaps.cmp_b2r ,transform=ccrs.PlateCarree(),extend = 'both') #RdYlBu_r nice colormap
    
    # contours = ax.contour(pco2_ml.LON, pco2_ml.LAT, pco2sen_ml[i],
    #                       colors='black', linewidths=0.6, levels= 29)
    # ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours
    ax.coastlines(zorder =11)
    ax.add_feature(cf.LAND, facecolor="lightgrey", zorder=10)
  
 
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

# Add vertical colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', label='\u00b5atm')
cbar.set_ticks= levels   # ticks every 10 units (integers)


plt.suptitle('INCOIS (ML) Seasonal climatology of pCO₂ in Indian Ocean (1980–2019) ', fontsize=16, y= 0.92)
plt.show()

#%%

