#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 00:37:43 2025

@author: bobco-08
"""

import numpy as np
import xarray as xr

file  = '/home/bobco-08/24cl05012/CO2/data/data_1/monthly_climatology_pco2_takahashi_2019'

clim_data  = np.genfromtxt(file)
clim_data = clim_data[1::]

#%%

lat_min,lat_max = -30,30
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
output_path = "/home/bobco-08/24cl05012/CO2/data/data_1/takahashi_clim_M_2019"
pco2_clim.to_netcdf(output_path)

print("Subset saved to:", output_path)

#%%
import numpy as np
import xarray as xr
clim_file = "/home/bobco-08/24cl05012/CO2/data/data_1/takahashi_clim_M_2019"

ds = xr.open_dataset(clim_file)

la = ds['lat']
lo = ds['lon']
mo = ds['month']
psw = ds['pco2_sw']

#%%

seasons = xr.full_like(mo, "", dtype=object)


seasons = xr.where((mo >= 3) & (mo <= 5), "MAM",seasons)
seasons = xr.where((mo >= 6) & (mo <= 9), "JJAS", seasons)
seasons = xr.where((mo >= 10) & (mo <= 11), "ON", seasons)
seasons = xr.where((mo == 12)|(mo <= 2), "DJF", seasons)

psw  = psw.assign_coords(custom_season=("month", seasons.data))

psw_season_mean = psw.groupby('custom_season').mean()

season_order = ["DJF", "MAM", "JJAS", "ON"]
psw_season_mean = psw_season_mean.reindex({"custom_season": season_order})

#%%
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean.cm as cm
import cmaps


fig, axs = plt.subplots(2, 2, figsize=(18,14),dpi = 600,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace':0.06,'hspace': -0.01})
season = ['DJF','MAM','JJAS',"ON"]
xticks = np.arange(35,120,10)
yticks = np.arange(-30,40,10)


season_lv = np.arange(300,460,5)
slv = np.arange(300,460,10)
#plotting

for i, ax in enumerate(axs.flat):
    im = ax.contourf(psw_season_mean.lon,psw_season_mean.lat,psw_season_mean[i], levels =season_lv,
                     cmap= cmaps.WhiteBlueGreenYellowRed,
                     transform=ccrs.PlateCarree(),extend = 'both') 
    ax.coastlines(zorder =15)
    ax.coastlines(zorder = 12)
    land = cf.NaturalEarthFeature(
        'physical', 'land', '10m',
        edgecolor='black',
        facecolor='gray'
    )
    ax.add_feature(land, zorder=2)

# subplot individual title

    ax.text(0.75, 0.98, f"Mean ({season[i]})",transform=ax.transAxes,         
        fontsize=18, zorder= 18,fontweight='bold',
        va='top', ha='right',)       

#contours  
 
    # contours = ax.contour(lon,lat,co2f_season[i],levels = season_lv,
    #                           colors='black', linewidths=0.6)
    # ax.clabel(contours, inline=True, fontsize=12, fmt="%.0f")        

# ticks and gridlines
# ticks and gridlines

    ax.set_xticks(xticks, crs=ccrs.PlateCarree()) # for adding ticks 
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    
    ax.tick_params(axis='both',
                   direction='out',   # makes arrow-like ticks
                   length=5,
                   width=1.2,
                   labelsize=20)
    
      # Hide labels depending on subplot index
    if i in [0, 1]:          # top row
        ax.tick_params(labelbottom=False)
    
    if i in [1, 3]:          # right column
        ax.tick_params(labelleft=False)

#colorbar

cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.70]) 
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical',ticks =slv )
cbar.ax.tick_params(axis='both', which='major',
               direction='out',   # ticks pointing outward
               length=6,          # tick length
               width=1,labelsize=20)


#titles
cbar.set_label('\u00b5atm',fontsize = 20,labelpad=12,
               fontweight = 'bold')
# plt.suptitle('Seasonal climatology of pCO2 in indian ocean 4° × 5° (1970 - 2007)', fontsize=20, y= 0.90)
# plt.show()
