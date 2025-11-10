#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 16:28:55 2025

@author: bobco-08
"""

import xarray as xr
import xesmf as xe
import numpy as np



# your dataset
ds = xr.open_dataset("/home/bobco-08/Desktop/24cl05012/CO2/data/ibr_pco2_mon_1980_2019.nc")

lat_res = ds['LAT']
lon_res = ds['LON']
time_res = ds['TIME']
pco2_res = ds['pCO2_Original']

pco2_IO_res = pco2_res.sel(LAT = slice(0,30),LON = slice(38,110))

# Save to new NetCDF
output_path = "/home/bobco-08/Desktop/24cl05012/CO2/data/pco2_IO_res.nc"
pco2_IO_res.to_netcdf(output_path)

print("Subset saved to:", output_path)
#%%

data = xr.open_dataset("/home/bobco-08/Desktop/24cl05012/CO2/data/pco2_IO_res.nc")


# build target grid
lat_1 = np.arange(1,30, 2)     # exactly the grid you showed
lon_1 = np.arange(38.75,111.3, 2.5)      # or whatever spacing you need
target_grid = xr.Dataset({"lat": (["lat"], lat_1),"lon": (["lon"], lon_1)})



# create a regridder and apply (bilinear, conservative, etc.)
regridder = xe.Regridder(data, target_grid, method="bilinear")


pco2_res_original  = data['pCO2_Original']

ocean_mask = xr.where(np.isfinite(pco2_res_original ), 1.0, 0.0)

numerator = regridder(pco2_res_original* ocean_mask)
denominator = regridder(ocean_mask)

pco2_regridded_coastal = (numerator / denominator).where(denominator > 0)
pco2_regridded_coastal.name = 'pCO2_regridded'

pco2_res_coastal = pco2_regridded_coastal

#%%
months = pco2_res_coastal['TIME'].dt.month

seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8]), "JJA", seasons)
seasons = xr.where(months.isin([9, 10, 11]), "SON", seasons)

pco2_io_original = pco2_res_coastal.assign_coords(custom_season=("TIME", seasons.data))

pco2_io_original_mean = pco2_io_original.groupby('custom_season').mean()

season_order = ["DJF", "MAM", "JJA", "SON"]
pco2_io_original_mean = pco2_io_original_mean.reindex({"custom_season": season_order})

#%%
a,b = np.nanpercentile(pco2_io_original_mean,[1,99])
print(a,b)
#%%

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf




fig, axs = plt.subplots(2, 2, figsize=(15,7),dpi = 200 ,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':0.03,'hspace': 0.010})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    levels = np.linspace(250,500,60)
    im = ax.contourf(pco2_io_original_mean.lon,pco2_io_original_mean.lat,pco2_io_original_mean[i], 
                     levels,cmap= 'seismic' ,transform=ccrs.PlateCarree()) #RdYlBu_r nice colormap
    
    contours = ax.contour(pco2_io_original_mean.lon,pco2_io_original_mean.lat,pco2_io_original_mean[i],
                          colors='black', linewidths=0.6, levels= 15)
    
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
        

# Add vertical colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', label='\u00b5atm')
cbar.set_ticks(np.arange(250,500, 25))   # ticks every 10 units (integers)
cbar.set_ticklabels([str(i) for i in range(250,500,25)])

plt.suptitle('Coarse resolution seasonal climatology of IBR model (1980-2019)', fontsize=16, y= 0.95)
plt.show()
