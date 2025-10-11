#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 00:37:43 2025

@author: bobco-08
"""

import numpy as np
import xarray as xr

file  = '/home/bobco-08/24cl05012/CO2/data/monthly_climatology_pco2'

clim_data  = np.genfromtxt(file)
clim_data = clim_data[1::]

#%%

lat_min,lat_max = -40,32
month_min,month_max = 1,12
lon_min,lon_max = 30,120


lat = clim_data[0:,0]
lon = clim_data[0:,1]
month = clim_data[0:,2]

mask = (lat >= lat_min) & (lat <= lat_max) & (month >= month_min) & (month <= month_max) & (lon >= lon_min) & (lon <= lon_max)

subset = clim_data[mask, :]

clim_io = subset[:,0:4]

uni_lat = np.unique_values(clim_io[0:,0])
uni_lon = np.unique_values(clim_io[0:,1])
uni_mon = np.unique_values(clim_io[0:,2])

# --- Prepare an empty 3D array ---

clim_grid = np.full((len(uni_mon), len(uni_lat), len(uni_lon)), np.nan)

# --- Fill data into the grid ---
for row in clim_io:
    la, lo, mo, val = row
    i = np.where(uni_mon == mo)[0][0]
    j = np.where(uni_lat == la)[0][0]
    k = np.where(uni_lon == lo)[0][0]
    clim_grid[i, j, k] = val
    
pco2_clim  = xr.Dataset(
    {
        "pco2_sw": (("month", "lat", "lon"), clim_grid)
    },
    coords={
        "lat": uni_lat,
        "lon": uni_lon,
        "month": uni_mon
    }
)    
#%%
output_path = "/home/bobco-08/24cl05012/CO2/data/pco2_climatology_taka.nc"
pco2_clim.to_netcdf(output_path)

print("Subset saved to:", output_path)

#%%
import numpy as np
import xarray as xr
clim_file = "/home/bobco-08/24cl05012/CO2/data/pco2_climatology_taka.nc"

ds = xr.open_dataset(clim_file)

la = ds['lat']
lo = ds['lon']
mo = ds['month']
psw = ds['pco2_sw']

#%%

seasons = xr.full_like(mo, "", dtype=object)


seasons = xr.where((mo >= 3) & (mo <= 5), "MAM",seasons)
seasons = xr.where((mo >= 6) & (mo <= 8), "JJA", seasons)
seasons = xr.where((mo >= 9) & (mo <= 11), "SON", seasons)
seasons = xr.where((mo == 12)|(mo <= 2), "DJF", seasons)

psw  = psw.assign_coords(custom_season=("month", seasons.data))

psw_season_mean = psw.groupby('custom_season').mean()

season_order = ["DJF", "MAM", "JJA", "SON"]
psw_season_mean = psw_season_mean.reindex({"custom_season": season_order})

#%%
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean.cm as cm


fig, axs = plt.subplots(2, 2, figsize=(15,10),dpi = 300 ,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':-0.18,'hspace': 0.06})
season = ['DJF','MAM','JJA',"SON"]


for i, ax in enumerate(axs.flat):
    levels = np.linspace(250,500,60)
    im = ax.contourf(psw_season_mean.lon,psw_season_mean.lat,psw_season_mean[i], 
                     levels,cmap= cm.balance ,transform=ccrs.PlateCarree()) #RdYlBu_r nice colormap
    
    contours = ax.contour(psw_season_mean.lon,psw_season_mean.lat,psw_season_mean[i],
                          colors='black', linewidths=0.6, levels= 10)
    
    ax.clabel(contours, inline=True, fontsize= 7, fmt="%.1f",levels=[contours.levels[0], contours.levels[2], contours.levels[3],contours.levels[5],contours.levels[6],contours.levels[7],contours.levels[9]])  # label the contours
    ax.coastlines(zorder =11)
    ax.add_feature(cf.LAND, facecolor="lightgrey", zorder=10)
  
 
# Used to add the subplot title inside the plot    
    ax.text(0.24, 0.97, f"{season[i]}",          # x, y in axes fraction (0–1)
        transform=ax.transAxes,         # interpret coords relative to axes
        fontsize=15, zorder= 18,fontweight='bold',fontfamily = 'serif',
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
cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', label='\u00b5atm')
cbar.set_ticks(np.arange(250,500, 25))   # ticks every 10 units (integers)
cbar.set_ticklabels([str(i) for i in range(250,500,25)])

plt.suptitle('Seasonal climatology of pCO2 in indian ocean  (4° × 5°)', fontsize=16, y= 0.92)
plt.show()
