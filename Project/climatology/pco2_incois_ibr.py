#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 18:54:13 2025

@author: bobco-08
"""

# Data loading

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean.cm as cm

file = "/home/bobco-08/24cl05012/CO2/data/ibr_pco2_mon_1980_2019.nc"

data = xr.open_dataset(file)

lat =  data['LAT']
lon = data['LON']
time = data['TIME']
pco2_ibr = data['pCO2_Original']

#%%%

#Trying to understand how the data is placed most of the values are arranged as NaN.

valid_counts = []
nan_counts = []

for t in range(len(pco2_ibr.TIME)):
    arr = pco2_ibr.isel(TIME=t).values
    valid = np.count_nonzero(~np.isnan(arr))
    nans = np.count_nonzero(np.isnan(arr))
    valid_counts.append(valid)
    nan_counts.append(nans)

valid_counts = np.array(valid_counts)
nan_counts = np.array(nan_counts)

# Print summary
print("PCO2 Data Validation Summary")
print("============================")
print(f"Total TIME steps (months): {len(pco2_ibr.TIME)}")
print(f"Min valid points in a month: {valid_counts.min()}")
print(f"Max valid points in a month: {valid_counts.max()}")
print(f"Mean valid points per month: {valid_counts.mean():.0f}")
print(f"Min NaN points in a month: {nan_counts.min()}")
print(f"Max NaN points in a month: {nan_counts.max()}")
print(f"Mean NaN points per month: {nan_counts.mean():.0f}")
print(f"Number of months with at least 1 valid point: {(valid_counts>0).sum()}")
print(f"Number of months with all points NaN: {(valid_counts==0).sum()}")
#%%
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

file = "/home/bobco-08/Downloads/6D32B6C6026814DD3A4D7CCF5CA31201_ferret_listing.nc"
data = xr.open_dataset(file)

pco2_ibr = data['pCO2_Original']

# --- 1️⃣ Select the time slice 1980–2020 ---
pco2_sel = pco2_ibr.sel(TIME=slice("1980-01-01", "2020-12-31"))

# --- 2️⃣ For each lat/lon, check if *all* time steps are finite (no NaN) ---
has_full_data = pco2_sel.notnull().all(dim="TIME")
# has_full_data is a 2-D (LAT,LON) boolean DataArray

# --- 3️⃣ Plot: red = complete data, white = missing at least once ---
fig = plt.figure(figsize=(10,6))
ax  = plt.axes(projection=ccrs.PlateCarree())

# Convert to integers so we can color (1 = full data, 0 = not)
mask_int = has_full_data.astype(int)

im = ax.pcolormesh(
    pco2_ibr['LON'],
    pco2_ibr['LAT'],
    mask_int,
    cmap=plt.cm.get_cmap('Reds', 2),   # 2 discrete colors
    vmin=0, vmax=1
)

# Add coastlines, etc.
ax.coastlines()
ax.add_feature(cf.LAND, facecolor="lightgrey")
ax.set_title("Grid cells with complete monthly pCO₂ data (1980–2020)")

# Custom colorbar: 0 = missing, 1 = full data
cbar = plt.colorbar(im, ax=ax, ticks=[0,1], shrink=0.7)
cbar.ax.set_yticklabels(['No full record','Full record'])

plt.show()


#%%
 
pco2_ibr = pco2_ibr.sel(LAT = slice(-30,30),LON = slice(32.5,110))

months = pco2_ibr['TIME'].dt.month

# Build season labels as DataArray
seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8]), "JJA", seasons)
seasons = xr.where(months.isin([9, 10, 11]), "SON", seasons)

# Assign as coordinate (must use .data since 'seasons' is a DataArray)
pco2_flux = pco2_ibr.assign_coords(
    custom_season=("TIME", seasons.data))

pco2sen_mean = pco2_flux.groupby('custom_season').mean()

# define the climatological order you want
season_order = ["DJF", "MAM", "JJA", "SON"]

pco2sen_mean = pco2sen_mean.reindex({"custom_season": season_order})



#%%

# Figuring out the outliers in the data
p5, p95 = np.nanpercentile(pco2sen_mean.values, [0,100])
print(p5, p95)
#%%%
fig, axs = plt.subplots(2, 2, figsize=(15,10),dpi = 300 ,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':-0.26,'hspace': 0.06})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    levels = np.linspace(250,500,60)
    im = ax.contourf(pco2_flux.LON, pco2_flux.LAT, pco2sen_mean[i], 
                     levels,cmap= cm.balance ,transform=ccrs.PlateCarree()) #RdYlBu_r nice colormap
    
    contours = ax.contour(pco2_flux.LON, pco2_flux.LAT, pco2sen_mean[i],
                          colors='black', linewidths=0.6, levels= 29)
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
cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', label='\u00b5atm')
cbar.set_ticks(np.arange(250,500, 25))   # ticks every 10 units (integers)
cbar.set_ticklabels([str(i) for i in range(250,500,25)])
cbar.set_ticklabels([str(i) for i in range(250,500,25)])



plt.suptitle('IBR model Seasonal climatology of pCO₂ in Indian Ocean (1980–2019) ', fontsize=16, y= 0.92)
plt.show()
#%%
pco2_flux_years = pco2_flux.assign_coords(year=("TIME", pco2_flux.TIME.dt.year.data))

pco2_season_year = pco2_flux_years.groupby(["year","custom_season"]).mean("TIME")

pco2_season_mean = pco2_season_year.mean("year")

pco2_season_mean = pco2_season_mean.transpose("custom_season", "LAT", "LON")

pco2_season_mean = pco2_season_mean.where(pco2_season_mean != 0)

pco2_season_mean = pco2_season_mean.reindex({"custom_season": season_order})

#%%

pco2_season_std = pco2_season_year.std("year")

pco2_season_std = pco2_season_std.transpose("custom_season", "LAT", "LON")

pco2_season_std = pco2_season_std.where(pco2_season_mean != 0)

pco2_season_std = pco2_season_std.reindex({"custom_season": season_order})


#%%

s5, s95 = np.nanpercentile(pco2_season_std.values, [1,99])
print(s5, s95)

#%%
fig, axs = plt.subplots(2, 2, figsize=(15,13),dpi = 150 ,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':0.05,'hspace': 0.05})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    levels = np.linspace(12,30,13)
    im = ax.contourf(pco2_flux.LON, pco2_flux.LAT, pco2_season_std[i], 
                     levels,cmap= 'viridis' ,transform=ccrs.PlateCarree()) #RdYlBu_r nice colormap
    
    contours = ax.contour(pco2_flux.LON, pco2_flux.LAT, pco2_season_std[i],
                          colors='black', linewidths=0.6, levels=20)
    ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours
    ax.coastlines(zorder =11)
    ax.add_feature(cf.LAND, facecolor="lightgrey", zorder=10)
  
 
# Used to add the subplot title inside the plot    
    ax.text(0.75, 0.98, f" Std ({season[i]})",          # x, y in axes fraction (0–1)
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



plt.suptitle('IBR model Standard deviation of Seasonal pCO₂ in Indian Ocean (1980-2019)', fontsize=16, y= 0.95)
#%%

fig, axs = plt.subplots(2, 4, figsize=(22, 10),dpi = 170, 
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace': 0.01, 'hspace': 0.0001})

season = ['DJF', 'MAM', 'JJA', "SON"]
levels = np.linspace(280,420,15)
# --- First row: Seasonal climatology ---'rainbow'
for i, ax in enumerate(axs[0, :]):
    im1 = ax.contourf(pco2_flux.LON, pco2_flux.LAT, pco2sen_mean[i], 
                      levels, cmap='viridis', transform=ccrs.PlateCarree(), extend='both')
    contours = ax.contour(pco2_flux.LON, pco2_flux.LAT, pco2sen_mean[i],
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
    im2 = ax.contourf(pco2_flux.LON, pco2_flux.LAT, pco2_season_std[i], 
                      levels=20, cmap='viridis', transform=ccrs.PlateCarree(), extend='both')
    contours = ax.contour(pco2_flux.LON, pco2_flux.LAT, pco2_season_std[i],
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
plt.suptitle('Seasonal climatology (top) and Standard deviation (bottom) of pCO₂  in Indian Ocean (1957–2023)',
             fontsize=16, y=0.95)

plt.show()

#%%

pco2 = pco2_ibr.sel(TIME = slice('1980-01-01','1980-12-31'))
