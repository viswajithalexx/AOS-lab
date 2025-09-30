#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 11:30:56 2025

@author: bobco-08
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 12:58:51 2025

@author: bobco-08
"""


import xesmf as xe
import numpy as np
import xarray as xr

# 1. LOAD DATA: Open the dataset and select the variable of interest
ds = xr.open_dataset("/home/bobco-08/Desktop/24cl05012/CO2/data/ibr_pco2_mon_1980_2019.nc")
pco2_original = ds['pCO2_Original']

# 2. DEFINE TARGET GRID: Create the new, coarser grid
lat_1 = np.arange(-30, 28, 2)
lon_1 = np.arange(30, 122.5, 2.5)
target_grid = xr.Dataset(
    {
        "lat": (["lat"], lat_1),
        "lon": (["lon"], lon_1),
    }
)

# 3. CREATE AN OCEAN MASK
# Create a mask where ocean cells are 1.0 and land cells (NaN) are 0.0
ocean_mask = xr.where(np.isfinite(pco2_original), 1.0, 0.0)

# 4. CREATE THE REGRIDDER
# We only need to create the regridder object once.
# Using a filename makes subsequent runs much faster.
regridder = xe.Regridder(
    ds,
    target_grid,
    method='conservative',
    filename='/home/bobco-08/Desktop/24cl05012/CO2/codes/conservative_weights.nc',
    reuse_weights=True
)

# 5. PERFORM MASKED REGRIDDING (THE KEY STEP)

# Numerator: Regrid the data multiplied by the mask (land values become zero)
numerator = regridder(pco2_original * ocean_mask)

# Denominator: Regrid the mask itself. This gives the fraction of each
# target cell that is covered by ocean in the source grid.
denominator = regridder(ocean_mask)

# Calculate the mean value for the ocean-only portion of each cell
# We use .where() to avoid division by zero for all-land cells.
pco2_regridded_coastal = (numerator / denominator).where(denominator > 0)
pco2_regridded_coastal.name = 'pCO2_regridded'


#%%
pco2_res = ds_coarse['pCO2_Original']

# Arranging months and assigning it as a custom season

months = pco2_res['TIME'].dt.month

# Build season labels as DataArray
seasons = xr.where(months.isin([12, 1, 2]), "DJF", "")
seasons = xr.where(months.isin([3, 4, 5]), "MAM", seasons)
seasons = xr.where(months.isin([6, 7, 8]), "JJA", seasons)
seasons = xr.where(months.isin([9, 10, 11]), "SON", seasons)

# Assign as coordinate (must use .data since 'seasons' is a DataArray)
pco2_flux_res = pco2_res.assign_coords(
    custom_season=("TIME", seasons.data))

pco2sen_mean_res = pco2_flux_res.groupby('custom_season').mean()

# define the climatological order you want
season_order = ["DJF", "MAM", "JJA", "SON"]

pco2sen_mean_res = pco2sen_mean_res.reindex({"custom_season": season_order})

#%%
#%%
fig, axs = plt.subplots(2, 2, figsize=(15,13),dpi = 150 ,subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace':0.05,'hspace': 0.05})
season = ['DJF','MAM','JJA',"SON"]

for i, ax in enumerate(axs.flat):
    levels = np.linspace(300,450,16)
    im = ax.contourf(pco2_flux_res.lon, pco2_flux_res.lat, pco2sen_mean_res[i], 
                     levels,cmap= 'viridis' ,transform=ccrs.PlateCarree()) #RdYlBu_r nice colormap
    
    contours = ax.contour(pco2_flux_res.lon, pco2_flux_res.lat, pco2sen_mean_res[i],
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
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', label='\u00b5atm')



plt.suptitle('IBR model Seasonal climatology of pCO₂ in Indian Ocean (1980–2019) ', fontsize=16, y= 0.95)
plt.show()
