#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 22:03:56 2025

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmaps


file = "/home/bobco-08/24cl05012/CO2/data/oc_v2025.flux.nc"

data = xr.open_dataset(file)

lat =  data['lat']
lon = data['lon']
time = data['mtime']
co2f = data['co2flux_ocean']*10**15
area = data['dxyp']

co2_flux = co2f/area #to convert per grid cell to per m^2

co2_io = co2_flux.sel(lat = slice(-30,30),lon = slice(30,120),mtime = slice('1960-01-01','2024-12-31'))
#%%%
months = co2_io['mtime'].dt.month

# Build season labels as DataArray
seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8, 9 ]), "JJAS", seasons)
seasons = xr.where(months.isin([10, 11]), "ON", seasons)


# Assign as coordinate (must use .data since 'seasons' is a DataArray)
co2_io = co2_io.assign_coords(custom_season=("mtime", seasons.data))

co2io_smean = co2_io.groupby('custom_season').mean()

co2io_smean = co2io_smean.where(co2io_smean != 0)

# define the climatological order you want
season_order = ["DJF", "MAM", "JJAS", "ON"]

co2io_smean = co2io_smean.reindex({"custom_season": season_order})
co2io_smeanv = co2io_smean.values

#%% Seasonal climatology of pCO₂ in Indian Ocean (1980–2019)
import matplotlib.ticker as mticker
fig, axs = plt.subplots(2, 2, figsize=(15,10),dpi = 300 ,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':0.02,'hspace': 0.06})
season = ['DJF','MAM','JJAS',"ON"]


for i, ax in enumerate(axs.flat):
    plot_levels = np.arange(-40,40.4,4)
    im = ax.contourf(co2_io.lon, co2_io.lat, co2io_smean[i], 
                     levels = plot_levels,cmap= 'RdBu_r' ,transform=ccrs.PlateCarree(),extend ='both' ) #RdYlBu_r nice colormap
    
    # contours = ax.contour(pco2_flux.lon, pco2_flux.lat, pco2sen_mean[i],
    #                       colors='black', linewidths=0.6, levels= 20)
    # ax.clabel(contours, inline=True, fontsize=8, fmt="%.0f")  # label the contours
    ax.coastlines(zorder =15)
    ax.add_feature(cf.LAND, facecolor="lightgrey", zorder=10)
  
 
# Used to add the subplot title inside the plot    
    ax.text(0.78, 0.98, f"Mean ({season[i]})",          # x, y in axes fraction (0–1)
        transform=ax.transAxes,         # interpret coords relative to axes
        fontsize=20,zorder= 18,fontweight='bold',
        va='top', ha='right',)            # vertical & horizontal alignment
    gl = ax.gridlines(draw_labels = True,linewidth = 0.5 , color = 'grey', alpha =0.1)
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.bottom_labels = False
    ax.tick_params(axis='both', which='major', 
                   direction='in',      # <--- This creates the "inward arrow" look
                   length=10, width=1, 
                   labelsize=15, 
                   top=False, right=False)  
    # 2. Logic for Outer Labels Only
    # If in Left Column (0 or 2), show Y-labels (Left)
    if i % 2 == 0: 
        gl.left_labels = True
    
    # If in Bottom Row (2 or 3), show X-labels (Bottom)
    if i >= 2: 
        gl.bottom_labels = True
    gl.xlabel_style = {'fontsize': 22, 'color': 'black'}
    gl.ylabel_style = {'fontsize': 22, 'color': 'black'}
  
# Add vertical colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(im, cax=cbar_ax ,ticks = plot_levels,orientation='vertical')
cbar.set_label(label='gC/$m^2$yr',fontsize= 20)
# cbar.set_ticks = plot_levels 
cbar.ax.tick_params(labelsize = 20)  
# plt.suptitle('Seasonal climatology of sea-air CO$_2$ flux in IO (1960–2024)', fontsize=21, y=0.95)
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/co2_flux/mean_co2flux_IO.pdf', dpi=600, bbox_inches='tight') 
plt.show()

#%%
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter # <--- Needed for °E / °N style

# ... (Your data loading code here) ...

fig, axs = plt.subplots(2, 2, figsize=(15, 10), dpi=300, 
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace':0.02,'hspace': 0.05})
season = ['DJF', 'MAM', 'JJAS', "ON"]
plot_levels = np.arange(-40, 40.4, 4)

for i, ax in enumerate(axs.flat):
    # 1. Plot Data
    im = ax.contourf(co2_io.lon, co2_io.lat, co2io_smean[i], 
                     levels=plot_levels, cmap='RdBu_r', 
                     transform=ccrs.PlateCarree(), extend='both')
    
    ax.coastlines(zorder=15)
    ax.add_feature(cf.LAND, facecolor="lightgrey", zorder=10)
    
    # 2. Set Extent (CRITICAL: Ensures your ticks are actually visible on screen)
    # Matches your requested ticks [40-100] and [-20 to 20]
    ax.set_extent([30, 120, -29, 30], crs=ccrs.PlateCarree())

    # 3. Add Title
    ax.text(0.78, 0.97, f"Mean ({season[i]})", transform=ax.transAxes, 
            fontsize=20, fontweight='bold', va='top', ha='right', zorder=18)

    # ---------------------------------------------------------
    # THE FIX: Use Standard Ticks instead of Gridlines Labels
    # ---------------------------------------------------------
    
    # A. Set the exact location of the ticks
    ax.set_xticks([40, 60, 80, 100], crs=ccrs.PlateCarree())
    ax.set_yticks([-20, -10, 0, 10, 20], crs=ccrs.PlateCarree())
    
    # B. Format them as Degrees (adds °E and °N)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    
    # C. Create the "Inward Arrows" (Ticks)
    ax.tick_params(axis='both', which='major', labelsize=15,
                   direction='out',      # <--- Points ticks INWARD
                   length=5, width=1.5, # <--- Size of the tick
                   top=False, right=False) # Turn off ticks on top/right borders
    
    # D. Handle Outer Labels Logic (Hide inner labels)
    if i % 2 != 0: # Right column (1, 3)
        ax.tick_params(labelleft=False) # Hide Y labels
        
    if i < 2:      # Top row (0, 1)
        ax.tick_params(labelbottom=False) # Hide X labels

    # ---------------------------------------------------------
    # OPTIONAL: Keep Gridlines ONLY for the faint inner grid
    # ---------------------------------------------------------
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, 
                      linewidth=0.5, color='grey', alpha=0.1, linestyle='--')


# --- Colorbar ---
cbar_ax = fig.add_axes([0.92, 0.14, 0.02, 0.7]) 
cbar = fig.colorbar(im, cax=cbar_ax, ticks=plot_levels[::2], orientation='vertical')
cbar.set_label(label='gC/$m^2$yr', fontsize=20)
cbar.ax.tick_params(labelsize=20)

# plt.suptitle('Seasonal climatology of sea-air CO$_2$ flux in Indian Ocean (1960–2024)', fontsize=2, y=0.94)
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/co2_flux/mean_co2flux_IO.pdf', dpi=600, bbox_inches='tight') 
plt.show()
#%%
a,b = np.nanpercentile(co2io_smean,([5,95]))



#%% for djf

co2flux_djf = co2_io.resample(mtime = 'QS-DEC').mean(dim=('mtime'))
co2flux_djf = co2flux_djf.isel(mtime=slice(1, None))
co2_flux_djf = co2flux_djf.isel(mtime = co2flux_djf.mtime.dt.month == 12)

co2_flux_djf = co2_flux_djf.where(co2_flux_djf != 0)

co2djf_sd = co2_flux_djf.std('mtime')
#%% for mam
co2flux_mam = co2_io.sel(mtime=co2_io.mtime.dt.month.isin([3,4,5]))
co2_flux_mam = co2flux_mam.resample(mtime='YE').mean(dim='mtime')
co2_flux_mam = co2_flux_mam.where(co2_flux_mam != 0)

co2mam_sd = co2_flux_mam.std('mtime')

#%% for jjas 
co2flux_jjas = co2_io.sel(mtime=co2_io.mtime.dt.month.isin([6,7,8,9]))
co2_flux_jjas = co2flux_mam.resample(mtime='YE').mean(dim='mtime')
co2_flux_jjas = co2_flux_mam.where(co2_flux_jjas != 0)

co2jjas_sd = co2_flux_jjas.std('mtime')

#%%
co2flux_so = co2_io.sel(mtime=co2_io.mtime.dt.month.isin([10,11]))
co2_flux_so = co2flux_so.resample(mtime='YE').mean(dim='mtime')
co2_flux_so = co2_flux_so.where(co2_flux_so != 0)

co2so_sd = co2_flux_so.std('mtime')

#%% Standard deviation of Seasonal pCO₂ in Indian Ocean (1980 -2019)

fig, axs = plt.subplots(2, 2, figsize=(15,10), dpi =300,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':0.02,'hspace': 0.06})
season = ['DJF','MAM','JJAS',"ON"]
season_sd = [co2djf_sd,co2mam_sd,co2jjas_sd,co2so_sd]

for i, ax in enumerate(axs.flat):
    sd_levels = np.arange(0,6.5,0.5)
    im = ax.contourf(co2_io.lon, co2_io.lat,season_sd[i], 
                  levels = sd_levels ,cmap= cmaps.WhiteBlueGreenYellowRed,transform=ccrs.PlateCarree(),extend ='max')
    # contours = ax.contour(pco2_flux.lon, pco2_flux.lat, pco2_season_std[i],
    #                       colors='black', linewidths=0.6, levels=15)
    # ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours

    ax.coastlines(zorder = 15)
    ax.add_feature(cf.LAND,color ='grey',zorder=11)
 
# Used to add the subplot title inside the plot    
    ax.text(0.75, 0.98, f"{season[i]} (SD)",          # x, y in axes fraction (0–1)
        transform=ax.transAxes,         # interpret coords relative to axes
        fontsize=12, zorder= 18,fontweight='bold',
        va='top', ha='right',)            # vertical & horizontal alignment
    gl = ax.gridlines(draw_labels = True,visible =False,linewidth = 0.5 , color = 'grey', alpha =0.5)
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'fontsize': 15}
    gl.ylabel_style = {'fontsize': 15}
    # remove *bottom labels* only for first row (i = 0,1)
    if i in [0, 1]:
        gl.bottom_labels = False
    if i in[1,3]:
        gl.left_labels = False
        

    

cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(im, cax=cbar_ax,ticks = sd_levels, orientation='vertical')
cbar.set_label(label='gC/$m^2$yr',fontsize = 20)
# cbar.set_ticks = levels_1   
cbar.ax.tick_params(labelsize = 15) 
plt.suptitle('Standard deviation of Seasonal sea-air CO$_2$ flux in IO (1960-2024)', fontsize=21, y=0.92)

plt.savefig('/home/bobco-08/24cl05012/CO2/plot/poster/co2_flux/co2_flux_sd', dpi=600, bbox_inches='tight') 
plt.show()

#%%

c,d = np.nanpercentile(co2so_sd,([5,]))
