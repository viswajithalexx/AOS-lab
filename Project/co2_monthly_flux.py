#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 17:35:35 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature
import cmaps



file = "/home/bobco-08/24cl05012/CO2/data/oc_v2025.pCO2 (1).nc"

data = xr.open_dataset(file,decode_timedelta=True)


lat =  data['lat']
lon = data['lon']
time = data['mtime']
co2_ocean = data['pCO2']

co2_flux = co2_ocean.sel(lat = slice(-30,30),lon = slice(30,110))


co2mon = co2_flux.groupby('mtime.month').mean(dim = 'mtime')
co2mon = co2mon.where(co2mon != 0)


#%%

fig, axs = plt.subplots(3, 4, figsize=(18, 14), subplot_kw={'projection': ccrs.PlateCarree()})

months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# 1. Your plot levels
plot_levels = np.arange(300,480, 10)

for i, ax in enumerate(axs.flat):
    
    # Use the new levels and add extend='both'
    im = ax.contourf(co2_flux.lon, co2_flux.lat, co2mon[i], 
                     levels=plot_levels, 
                     cmap=cmaps.cmp_b2r,  # CHANGED: Better colormap for sequential data
                     extend='both',       # ADDED: To handle values outside levels
                     transform=ccrs.PlateCarree())
    
    ax.text(0.57, 0.89, months[i], 
            transform=ax.transAxes,  # Use axes coordinates (0-1)
            fontsize=12, 
            fontweight='bold',
            color='black',
            zorder = 20,
            )
    ax.coastlines(zorder=15)
    # 2. Add land and borders for better context
    ax.add_feature(cfeature.LAND, color='grey', zorder=5)

    # 3. Use ax.gridlines for cleaner label control
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, visible=False)
    gl.top_labels = False
    gl.right_labels = False
    # Only show labels on the outer-most plots
    gl.left_labels = (i % 4 == 0) # Only first column
    gl.bottom_labels = (i >= 8)   # Only bottom row

cbar_ax = fig.add_axes([1.02, 0.30, 0.015, 0.6]) 

# 2. Draw the colorbar in that new axis
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label(label='$\mu$atm', labelpad=10)
# 4. Corrected main title to match your data
plt.suptitle('Monthly Mean pCO₂ in Indian Ocean (1957–2024)', fontsize=16)

# 5. FIXED: Changed bottom rect from 0.3 to 0.03 to remove the gap
plt.tight_layout(rect=[0, 0.2, 1.0, 1.01])
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/monthly_mean_pco2_IO.tiff', dpi=500, bbox_inches='tight') 
plt.show()