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

pco2_flux = pco2.sel(lat = slice(0,30,),lon = slice(38,80),mtime = slice('1980-01-01','2019-12-31'))


#%%
months = pco2_flux['mtime'].dt.month

# Build season labels as DataArray
seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8]), "JJA", seasons)
seasons = xr.where(months.isin([9, 10, 11]), "SON", seasons)


# Assign as coordinate (must use .data since 'seasons' is a DataArray)
pco2_flux = pco2_flux.assign_coords(custom_season=("mtime", seasons.data))

pco2sen_mean = pco2_flux.groupby('custom_season').mean()

pco2sen_mean = pco2sen_mean.where(pco2sen_mean != 0)

# define the climatological order you want
season_order = ["DJF", "MAM", "JJA", "SON"]

pco2sen_mean = pco2sen_mean.reindex({"custom_season": season_order})

#%%
p5,p90 = np.nanpercentile(pco2sen_mean.values,[5,95])
print(p5,p90)

#%% Seasonal climatology of pCO₂ in Indian Ocean (1980–2019)

fig, axs = plt.subplots(2, 2, figsize=(15,13),dpi = 150 ,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':0.15,'hspace': 0.19})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    levels = np.linspace(300,500,16)
    im = ax.contourf(pco2_flux.lon, pco2_flux.lat, pco2sen_mean[i], 
                     levels,cmap= 'rainbow' ,transform=ccrs.PlateCarree()) #RdYlBu_r nice colormap
    
    contours = ax.contour(pco2_flux.lon, pco2_flux.lat, pco2sen_mean[i],
                          colors='black', linewidths=0.6, levels= 20)
    ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours
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
        
# Add colorbar
cbar_ax = fig.add_axes([0.25, 0.50, 0.5, 0.015]) # left,bottom,width,height
cbar = fig.colorbar(im, cax=cbar_ax,orientation = 'horizontal', label= '\u00b5atm')

plt.suptitle('Seasonal climatology of pCO₂ in Indian Ocean (1980–2019)', fontsize=16, y= 0.95)
plt.show()

#%% Monthly Mean pCO₂ (1980–2019)

pco2mon = pco2_flux.groupby('mtime.month').mean(dim = 'mtime')
pco2mon = pco2mon.where(pco2mon != 0)



fig, axs = plt.subplots(3, 4, figsize=(25,24), subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':-0.55,'hspace': 0.02})
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

for i, ax in enumerate(axs.flat):
    
    im = ax.contourf(pco2_flux.lon, pco2_flux.lat, pco2mon[i], 
                     levels=20 , cmap='PiYG', transform=ccrs.PlateCarree())
    # contours = ax.contour(pco2_flux.lon, pco2_flux.lat, pco2mon[i],
    #                       colors='black', linewidths=0.6, levels=8)
    # ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours
    ax.coastlines()
    ax.add_feature(cf.BORDERS, linewidth=0.5)
    ax.add_feature(cf.LAND,linewidth=0.5, zorder =11)
    
    
    # Used to add the subplot title inside the plot    
    ax.text(0.75, 0.98, months[i],          # x, y in axes fraction (0–1)
            transform=ax.transAxes,         # interpret coords relative to axes
            fontsize=12, zorder= 18,fontweight='bold',
            va='top', ha='right',)            # vertical & horizontal alignment
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.5)
    gl.top_labels =  False
    gl.right_labels =  False
    
    if i in [0,1,2,3,4,5,6,7]:
        gl.bottom_labels = False
    if i in[1,2,3,5,6,7,9,10,11,12]:
        gl.left_labels = False
        
    
# Add colorbar
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.85, 0.25, 0.015, 0.5]) # left,bottom,width,height
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label(label =r'\u00b5atm',labelpad = 10)


plt.suptitle('Monthly Mean pCO₂ (1980–2019)', fontsize=16,y = 0.95)
plt.show()

#%%

pco2_flux_years = pco2_flux.assign_coords(year=("mtime", pco2_flux.mtime.dt.year.data))

pco2_season_year = pco2_flux_years.groupby(["year","custom_season"]).mean("mtime")

pco2_season_mean = pco2_season_year.mean("year")

pco2_season_mean = pco2_season_mean.transpose("custom_season", "lat", "lon")

pco2_season_mean = pco2_season_mean.where(pco2_season_mean != 0)

pco2_season_mean = pco2_season_mean.reindex({"custom_season": season_order})

#%%

pco2_season_std = pco2_season_year.std("year")

pco2_season_std = pco2_season_std.transpose("custom_season", "lat", "lon")

pco2_season_std = pco2_season_std.where(pco2_season_mean != 0)

pco2_season_std = pco2_season_std.reindex({"custom_season": season_order})



#%% Standard deviation of Seasonal pCO₂ in Indian Ocean (1980 -2019)


fig, axs = plt.subplots(2, 2, figsize=(15,13), dpi =150,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':-0.15,'hspace': 0.19})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    
    levels_1 = np.linspace(12,30,13)
    im = ax.contourf(pco2_flux.lon, pco2_flux.lat,pco2_season_std[i], 
                     levels_1,cmap= 'viridis',transform=ccrs.PlateCarree())
    contours = ax.contour(pco2_flux.lon, pco2_flux.lat, pco2_season_std[i],
                          colors='black', linewidths=0.6, levels=15)
    ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours

    ax.coastlines()
    ax.add_feature(cf.BORDERS, linewidth=0.5)
    ax.add_feature(cf.LAND,color ='grey',zorder=11)
 
# Used to add the subplot title inside the plot    
    ax.text(0.75, 0.98, f"std ({season[i]})",          # x, y in axes fraction (0–1)
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



plt.suptitle('Standard deviation of Seasonal pCO₂ in Indian Ocean (1980 -2019)', fontsize=16, y= 0.95)
plt.show()

#%% Seasonal climatology (top) and Standard deviation (bottom) of pCO₂  in Indian Ocean (1980 -2019)

fig, axs = plt.subplots(2, 4, figsize=(22, 10),dpi = 170, 
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace': 0.01, 'hspace': 0.0001})

season = ['DJF', 'MAM', 'JJA', "SON"]
levels = np.linspace(280,420,15)
# --- First row: Seasonal climatology ---'rainbow'
for i, ax in enumerate(axs[0, :]):
    im1 = ax.contourf(pco2_flux.lon, pco2_flux.lat, pco2sen_mean[i], 
                      levels, cmap='seismic', transform=ccrs.PlateCarree(), extend='both')
    contours = ax.contour(pco2_flux.lon, pco2_flux.lat, pco2sen_mean[i],
                          colors='black', linewidths=0.6, levels=15)
    ax.clabel(contours, inline=True, fontsize=7, fmt="%.1f")

    ax.coastlines(zorder=11)
    
    ax.add_feature(cf.LAND, color='lightgrey', zorder=10)

    # 👉 Just the season name
    ax.text(0.75, 0.98, f"Mean ({season[i]})", transform=ax.transAxes,
            fontsize=11,zorder = 19, fontweight='bold', va='top', ha='right')

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.5)
    gl.right_labels, gl.top_labels = False, False
    if i != 0: gl.left_labels = False
    gl.bottom_labels = False   # remove bottom labels for top row


# --- Second row: Standard deviation ---
for i, ax in enumerate(axs[1, :]):
    im2 = ax.contourf(pco2_flux.lon, pco2_flux.lat, pco2_season_std[i], 
                      levels=20, cmap='seismic', transform=ccrs.PlateCarree(), extend='both')
    contours = ax.contour(pco2_flux.lon, pco2_flux.lat, pco2_season_std[i],
                          colors='black', linewidths=0.6, levels=15)
    ax.clabel(contours, inline=True, fontsize=7, fmt="%.1f")

    ax.coastlines(zorder=11)
   
    ax.add_feature(cf.LAND, color='lightgrey', zorder=10)

    # 👉 "std_" + season name
    ax.text(0.75, 0.98, f"std ({season[i]})", transform=ax.transAxes,
            fontsize=12,zorder = 19, fontweight='bold', va='top', ha='right')

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.5)
    gl.right_labels, gl.top_labels = False, False
    if i != 0: gl.left_labels = False


# --- Colorbars on the right side ---
# Climatology (top row)
cbar_ax1 = fig.add_axes([0.92, 0.52, 0.015, 0.35])  # [left, bottom, width, height]
cbar1 = fig.colorbar(im1, cax=cbar_ax1, orientation='vertical', label='\u00b5atm')

# Standard deviation (bottom row)
cbar_ax2 = fig.add_axes([0.92, 0.11, 0.015, 0.35])
cbar2 = fig.colorbar(im2, cax=cbar_ax2, orientation='vertical', label='\u00b5atm')

# --- Titles ---
plt.suptitle('Seasonal climatology (top) and Standard deviation (bottom) of pCO₂  in Indian Ocean (1980 -2019)',
             fontsize=16, y=0.95)

plt.show()

