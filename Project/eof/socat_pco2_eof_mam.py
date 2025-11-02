#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 11:18:01 2025

@author: bobco-08
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 21:19:34 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
from eofs.xarray import Eof
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean.cm as cm
import cmaps

file1 = "/home/bobco-08/24cl05012/CO2/data/oc_v2025.pCO2.nc"
data1 = xr.open_dataset(file1)

lat1 = data1['lat']
lon1 = data1['lon']
time1 = data1['mtime']
pco2_1  = data1['pCO2']

pco2_1 = pco2_1.where(pco2_1 != 0) 
#%%
socat = pco2_1.sel(lat = slice(-30,30),lon = slice(30,110))
#%%
da = socat

da_djf = da.sel(mtime=da.mtime.dt.month.isin([12, 1, 2]))
da_mam = da.sel(mtime=da.mtime.dt.month.isin([3,4,5]))
da_jjas = da.sel(mtime=da.mtime.dt.month.isin([6, 7, 8, 9]))
da_on = da.sel(mtime=da.mtime.dt.month.isin([10,11]))


#%%
da_djf_yearly_mean = da_djf.resample(mtime='Y').mean(dim='mtime')

da_djf_yearly_mean = da_djf_yearly_mean.dropna(dim='mtime', how='all')
djf_climatology = da_djf_yearly_mean.mean(dim='mtime')
djf_ano = da_djf_yearly_mean - djf_climatology

#%%
import scipy.signal as sig
import time

djf_ano = djf_ano.fillna(0)
djf_detrend = sig.detrend(djf_ano, axis=0)

djf_detrend[djf_detrend == 0] = np.nan

djf_detrend = djf_detrend[2:-2,:,:]
#%%
import eofs
from eofs.standard import Eof

# Assuming 'djf_detrend' is your 3D NumPy array (time, lat, lon)
# and 'djf_ano' is the xarray DataArray from which 'lat' is taken.
# Get the shape from the detrended data
n_time, n_lat, n_lon = djf_detrend.shape
# --- 1. Calculate the weights ---
# Use the xarray object 'djf_ano' which still has the latitude coordinates
# The correct weight for EOFs is the square-root of the cosine of latitude
# Make sure to convert degrees to radians for np.cos
eof_weights_lat = np.sqrt(np.cos(np.deg2rad(djf_ano['lat'].values)))
# --- 2. Reshape your DATA for the solver ---
# The solver needs 2D (time, space). Space is lat * lon.
data_2d = djf_detrend.reshape(n_time, n_lat * n_lon)
# --- 3. Reshape your WEIGHTS for the solver ---
# Weights must be 1D (space), matching the 2nd dim of the data.
# a. Broadcast the 1D (lat) weights to a 2D (lat, lon) grid
weights_2d_grid = np.broadcast_to(eof_weights_lat[:, np.newaxis], (n_lat, n_lon))
# b. Flatten the 2D grid to a 1D (space) array
weights_1d_space = weights_2d_grid.flatten()
# --- 4. Call the solver with the corrected 2D data and 1D weights ---
# Note: The data being passed to the Eof class *must* be a NumPy array.
solver = Eof(data_2d, weights=weights_1d_space)
#%%
neofs = 3
covmaps = solver.eofsAsCovariance(neofs=neofs)
#%%
covmaps_mode1 = covmaps[0].reshape(n_lat,n_lon)
covmaps_mode2 = covmaps[1].reshape(n_lat,n_lon)
covmaps_mode3 = covmaps[2].reshape(n_lat,n_lon)
#%%
plt.figure(figsize =(15,10),dpi = 300)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im = plt.pcolormesh(da.lon, da.lat, covmaps_mode1,vmin = -8,vmax =8,cmap=cmaps.cmp_b2r)
ax.add_feature(cf.LAND, color='grey', zorder=12)
ax.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
cbar = plt.colorbar(im, orientation='vertical', shrink=0.7)
cbar.set_label('$\mu$atm',labelpad = 5)

# Setting the suggested title label
ax.set_title('IO:DJF:1959-2022', fontsize=16,pad =20)
# ➡️ Add the Subtitle for Variance Explained using ax.text() ⬅️
ax.text(
    1.09,                                         # Center position (0.5 is the center of the x-axis)
    0.92,                                        # Vertical position (0.95 is just below the top edge)
    '63.04 % ', # The text string
    horizontalalignment='center',                # Center the text
    transform=ax.transAxes,                      # Use axes coordinates (0 to 1)
    fontsize=12,
    color='k'
    )
ax.text(
    0.024,                                         # Center position (0.5 is the center of the x-axis)
    1.02,                                        # Vertical position (0.95 is just below the top edge)
    'EOF 1 ', # The text string
    horizontalalignment='center',                # Center the text
    transform=ax.transAxes,                      # Use axes coordinates (0 to 1)
    fontsize=12,
    color='k'
)
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/eof/DJF/IO_DJF_19592022_eof1.png', dpi=500, bbox_inches='tight')
plt.show()
#%%
a,b = np.nanpercentile(covmaps_mode1,[5,95])
#%%
eofvar = solver.varianceFraction(neigs=neofs) * 100
eofvar
#%%
# Assuming 'solver' is your eofs.standard.Eof object
neofs = 3
# --- 1. Get PC Time Series ---
# pcs has shape (time, neofs) by default, using .T makes it (neofs, time).
# We keep the shape (time, neofs) for easier plotting, so remove the .T
pcs = solver.pcs(pcscaling=1, npcs=neofs)
#%% DJF EOF1

years = np.arange(1959,2023,1)

# --- 3. Plot the First PC (Mode 1) ---
plt.figure(figsize=(15,6), dpi=500) 

# Use the 'years' array for the x-axis
plt.plot(years, pcs[:, 0], linewidth=2.5, color='darkblue', label='PC1')

# Optional: Add a horizontal zero line for reference
plt.axhline(0, color='grey', linestyle='--', linewidth=1.0) 
plt.fill_between(
    years, pcs[:,0], 0,
    where= pcs[:,0] >= 0, # Condition: fill where PC1 data is greater than or equal to 0
    facecolor='red',     # Fill color
    alpha=0.3,           # Transparency
    interpolate=True,    # Ensures smooth filling even if line crosses zero between points
    zorder=1             # Keep behind the line and zero axis
)
# ➡️ Add shading for values BELOW 0 (blue) ⬅️
plt.fill_between(
    years, pcs[:,0], 0,
    where=pcs[:,0] < 0,  # Condition: fill where PC1 data is less than 0
    facecolor='blue',    # Fill color
    alpha=0.3,           # Transparency
    interpolate=True,    # Ensures smooth filling
    zorder=1             # Keep behind the line and zero axis
)
plt.xlim(years.min(),years.max())
plt.ylim(-3,3)
plt.ylabel('Standardized', fontsize=12)
plt.xlabel('Year (DJF Season)', fontsize=12)
plt.title('EOF 1 63.04%', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/eof/DJF/IO_DJF_19592022_pcs1.pdf', dpi=500, bbox_inches='tight')
plt.show()
#%%

exact_yr = covmaps_mode1*pcs[6,0]

plt.figure()
plt.figure(figsize =(15,10),dpi = 300)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines(zorder = 15)

im = plt.contourf(da.lon, da.lat,exact_yr, cmap=cmaps.cmp_b2r, extend='both')
ax.add_feature(cf.LAND, color='grey', zorder=12)
ax.gridlines(visible=True, draw_labels=True, dms=True, x_inline=False, y_inline=False)
cbar = plt.colorbar(im, orientation='vertical', shrink=0.7)
cbar.set_label('EOF-1 Spatial Pattern ($\mu$atm)',labelpad =20)
#%% djf_eof2
plt.figure(figsize =(15,10),dpi = 300)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im = plt.pcolormesh(da.lon, da.lat, covmaps_mode2,vmin = -8,vmax =8,cmap=cmaps.cmp_b2r)
ax.add_feature(cf.LAND, color='grey', zorder=12)
ax.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
cbar = plt.colorbar(im, orientation='vertical', shrink=0.7)
cbar.set_label('$\mu$atm',labelpad = 5)

# Setting the suggested title label
ax.set_title('IO:DJF:1959-2022', fontsize=16,pad =20)
# ➡️ Add the Subtitle for Variance Explained using ax.text() ⬅️
ax.text(
    1.09,                                         # Center position (0.5 is the center of the x-axis)
    0.92,                                        # Vertical position (0.95 is just below the top edge)
    '7.72 % ', # The text string
    horizontalalignment='center',                # Center the text
    transform=ax.transAxes,                      # Use axes coordinates (0 to 1)
    fontsize=12,
    color='k'
    )
ax.text(
    0.024,                                         # Center position (0.5 is the center of the x-axis)
    1.02,                                        # Vertical position (0.95 is just below the top edge)
    'EOF 2 ', # The text string
    horizontalalignment='center',                # Center the text
    transform=ax.transAxes,                      # Use axes coordinates (0 to 1)
    fontsize=12,
    color='k'
)
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/eof/DJF/IO_DJF_19592022_eof2.png', dpi=500, bbox_inches='tight')
plt.show()
#%%
years = np.arange(1959,2023,1)

# --- 3. Plot the First PC (Mode 1) ---
plt.figure(figsize=(15,6), dpi=300) 

# Use the 'years' array for the x-axis
plt.plot(years, pcs[:, 1], linewidth=2.5, color='darkblue', label='PC2')

# Optional: Add a horizontal zero line for reference
plt.axhline(0, color='grey', linestyle='--', linewidth=1.0) 
plt.fill_between(
    years, pcs[:,1], 0,
    where= pcs[:,1] >= 0, # Condition: fill where PC1 data is greater than or equal to 0
    facecolor='red',     # Fill color
    alpha=0.3,           # Transparency
    interpolate=True,    # Ensures smooth filling even if line crosses zero between points
    zorder=1             # Keep behind the line and zero axis
)
# ➡️ Add shading for values BELOW 0 (blue) ⬅️
plt.fill_between(
    years, pcs[:,1], 0,
    where=pcs[:,1] < 0,  # Condition: fill where PC1 data is less than 0
    facecolor='blue',    # Fill color
    alpha=0.3,           # Transparency
    interpolate=True,    # Ensures smooth filling
    zorder=1             # Keep behind the line and zero axis
)
plt.xlim(years.min(),years.max())
plt.ylim(-3,3)
plt.ylabel('Standardized', fontsize=12)
plt.xlabel('Year (DJF Season)', fontsize=12)
plt.title('EOF 2 7.72%', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/eof/DJF/IO_DJF_19592022_pcs2.pdf', dpi=500, bbox_inches='tight')
plt.show()
#%%
plt.figure(figsize =(15,10),dpi = 300)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im = plt.pcolormesh(da.lon, da.lat, covmaps_mode3,vmin = -8,vmax =8,cmap=cmaps.cmp_b2r)
ax.add_feature(cf.LAND, color='grey', zorder=12)
ax.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
cbar = plt.colorbar(im, orientation='vertical', shrink=0.7)
cbar.set_label('$\mu$atm',labelpad = 5)

# Setting the suggested title label
ax.set_title('IO:DJF:1959-2022', fontsize=16,pad =20)
# ➡️ Add the Subtitle for Variance Explained using ax.text() ⬅️
ax.text(
    1.09,                                         # Center position (0.5 is the center of the x-axis)
    0.92,                                        # Vertical position (0.95 is just below the top edge)
    '4.09 % ', # The text string
    horizontalalignment='center',                # Center the text
    transform=ax.transAxes,                      # Use axes coordinates (0 to 1)
    fontsize=12,
    color='k'
    )
ax.text(
    0.024,                                         # Center position (0.5 is the center of the x-axis)
    1.02,                                        # Vertical position (0.95 is just below the top edge)
    'EOF 3 ', # The text string
    horizontalalignment='center',                # Center the text
    transform=ax.transAxes,                      # Use axes coordinates (0 to 1)
    fontsize=12,
    color='k'
)
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/eof/DJF/IO_DJF_19592022_eof3.png', dpi=500, bbox_inches='tight')
plt.show()
#%%
# --- 3. Plot the First PC (Mode 1) ---
plt.figure(figsize=(15, 6), dpi=300) 

# Use the 'years' array for the x-axis
plt.plot(years, pcs[:, 2], linewidth=2.5, color='darkblue', label='PC3')

# Optional: Add a horizontal zero line for reference
plt.axhline(0, color='grey', linestyle='--', linewidth=1.0) 
plt.fill_between(
    years, pcs[:,2], 0,
    where= pcs[:,2] >= 0, # Condition: fill where PC1 data is greater than or equal to 0
    facecolor='red',     # Fill color
    alpha=0.3,           # Transparency
    interpolate=True,    # Ensures smooth filling even if line crosses zero between points
    zorder=1             # Keep behind the line and zero axis
)
# ➡️ Add shading for values BELOW 0 (blue) ⬅️
plt.fill_between(
    years, pcs[:,2], 0,
    where=pcs[:,2] < 0,  # Condition: fill where PC1 data is less than 0
    facecolor='blue',    # Fill color
    alpha=0.3,           # Transparency
    interpolate=True,    # Ensures smooth filling
    zorder=1             # Keep behind the line and zero axis
)
plt.xlim(years.min(),years.max())
plt.ylim(-3,3)
plt.ylabel('Standardized', fontsize=12)
plt.xlabel('Year (DJF Season)', fontsize=12)
plt.title('EOF 3 4.09%', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/eof/DJF/IO_DJF_19592022_pcs3.pdf', dpi=500, bbox_inches='tight')
plt.show()
#%%
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmaps

# --- Assumed Variables ---
# Make sure the following variables are defined from your previous cells:
# da, covmaps_mode1, covmaps_mode2, covmaps_mode3, pcs
# -----------------------------------

# --- 1. Create the Main Figure ---
# A 3x2 grid, sized to hold all plots
# Increased width for 2 columns (15*2) and height for 3 rows (10*3)
fig = plt.figure(figsize=(30, 30), dpi=300)

# Define the years array once
years = np.arange(1959, 2023, 1)

# ---------------------------------
# --- Row 1: EOF 1 and PC 1 ---
# ---------------------------------

# --- AX 1: EOF 1 Map (Top-Left) ---
ax1 = fig.add_subplot(3, 2, 1, projection=ccrs.PlateCarree())
ax1.coastlines(zorder=15)
im1 = ax1.pcolormesh(da.lon, da.lat, covmaps_mode1, vmin=-8, vmax=8, cmap=cmaps.cmp_b2r, transform=ccrs.PlateCarree())
ax1.add_feature(cf.LAND, color='grey', zorder=12)
gl1 = ax1.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl1.top_labels = False    # Added to clean up subplot
gl1.right_labels = False  # Added to clean up subplot
cbar1 = plt.colorbar(im1, ax=ax1, orientation='vertical', shrink=0.7) # Use ax=ax1
cbar1.set_label('$\mu$atm', labelpad=5)
ax1.set_title('IO:DJF:1959-2022', fontsize=16, pad=20)
ax1.text(
    1.09, 0.90, '63.04 % ', # Your original text
    horizontalalignment='center',
    transform=ax1.transAxes,
    fontsize=12,
    color='k'
)
ax1.text(
    0.024, 1.02, 'EOF 1 ', # Your original text
    horizontalalignment='center',
    transform=ax1.transAxes,
    fontsize=12,
    color='k'
)

# --- AX 2: PC 1 Plot (Top-Right) ---
ax2 = fig.add_subplot(3, 2, 2)
ax2.plot(years, pcs[:, 0], linewidth=2.5, color='darkblue', label='PC1')
ax2.axhline(0, color='grey', linestyle='--', linewidth=1.0)
ax2.fill_between(
    years, pcs[:, 0], 0,
    where=pcs[:, 0] >= 0,
    facecolor='red', alpha=0.3, interpolate=True, zorder=1
)
ax2.fill_between(
    years, pcs[:, 0], 0,
    where=pcs[:, 0] < 0,
    facecolor='blue', alpha=0.3, interpolate=True, zorder=1
)
ax2.set_xlim(years.min(), years.max())
ax2.set_ylim(-3, 3)
ax2.set_ylabel('Standardized', fontsize=12)
ax2.set_xlabel('Year (DJF Season)', fontsize=12)
ax2.set_title('EOF 1 63.04%', fontsize=14)
ax2.grid(True, linestyle=':', alpha=0.7)
ax2.legend(loc='upper right')

# ---------------------------------
# --- Row 2: EOF 2 and PC 2 ---
# ---------------------------------

# --- AX 3: EOF 2 Map (Middle-Left) ---
ax3 = fig.add_subplot(3, 2, 3, projection=ccrs.PlateCarree())
ax3.coastlines(zorder=15)
im3 = ax3.pcolormesh(da.lon, da.lat, covmaps_mode2, vmin=-8, vmax=8, cmap=cmaps.cmp_b2r, transform=ccrs.PlateCarree())
ax3.add_feature(cf.LAND, color='grey', zorder=12)
gl3 = ax3.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl3.top_labels = False    # Added to clean up subplot
gl3.right_labels = False  # Added to clean up subplot
cbar3 = plt.colorbar(im3, ax=ax3, orientation='vertical', shrink=0.7) # Use ax=ax3
cbar3.set_label('$\mu$atm', labelpad=5)
ax3.set_title('IO:DJF:1959-2022', fontsize=16, pad=20)
ax3.text(
    1.09, 0.90, '7.72 % ', # Your original text
    horizontalalignment='center',
    transform=ax3.transAxes,
    fontsize=12,
    color='k'
)
ax3.text(
    0.024, 1.02, 'EOF 2 ', # Your original text
    horizontalalignment='center',
    transform=ax3.transAxes,
    fontsize=12,
    color='k'
)

# --- AX 4: PC 2 Plot (Middle-Right) ---
ax4 = fig.add_subplot(3, 2, 4)
ax4.plot(years, pcs[:, 1], linewidth=2.5, color='darkblue', label='PC2')
ax4.axhline(0, color='grey', linestyle='--', linewidth=1.0)
ax4.fill_between(
    years, pcs[:, 1], 0,
    where=pcs[:, 1] >= 0,
    facecolor='red', alpha=0.3, interpolate=True, zorder=1
)
ax4.fill_between(
    years, pcs[:, 1], 0,
    where=pcs[:, 1] < 0,
    facecolor='blue', alpha=0.3, interpolate=True, zorder=1
)
ax4.set_xlim(years.min(), years.max())
ax4.set_ylim(-3, 3)
ax4.set_ylabel('Standardized', fontsize=12)
ax4.set_xlabel('Year (DJF Season)', fontsize=12)
ax4.set_title('EOF 2 7.72%', fontsize=14)
ax4.grid(True, linestyle=':', alpha=0.7)
ax4.legend(loc='upper right')

# ---------------------------------
# --- Row 3: EOF 3 and PC 3 ---
# ---------------------------------

# --- AX 5: EOF 3 Map (Bottom-Left) ---
ax5 = fig.add_subplot(3, 2, 5, projection=ccrs.PlateCarree())
ax5.coastlines(zorder=15)
im5 = ax5.pcolormesh(da.lon, da.lat, covmaps_mode3, vmin=-8, vmax=8, cmap=cmaps.cmp_b2r, transform=ccrs.PlateCarree())
ax5.add_feature(cf.LAND, color='grey', zorder=12)
gl5 = ax5.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl5.top_labels = False    # Added to clean up subplot
gl5.right_labels = False  # Added to clean up subplot
cbar5 = plt.colorbar(im5, ax=ax5, orientation='vertical', shrink=0.7) # Use ax=ax5
cbar5.set_label('$\mu$atm', labelpad=5)
ax5.set_title('IO:DJF:1959-2022', fontsize=16, pad=20)
ax5.text(
    1.09, 0.90, '4.09 % ', # Your original text
    horizontalalignment='center',
    transform=ax5.transAxes,
    fontsize=12,
    color='k'
)
ax5.text(
    0.024, 1.02, 'EOF 3 ', # Your original text
    horizontalalignment='center',
    transform=ax5.transAxes,
    fontsize=12,
    color='k'
)

# --- AX 6: PC 3 Plot (Bottom-Right) ---
ax6 = fig.add_subplot(3, 2, 6)
ax6.plot(years, pcs[:, 2], linewidth=2.5, color='darkblue', label='PC3')
ax6.axhline(0, color='grey', linestyle='--', linewidth=1.0)
ax6.fill_between(
    years, pcs[:, 2], 0,
    where=pcs[:, 2] >= 0,
    facecolor='red', alpha=0.3, interpolate=True, zorder=1
)
ax6.fill_between(
    years, pcs[:, 2], 0,
    where=pcs[:, 2] < 0,
    facecolor='blue', alpha=0.3, interpolate=True, zorder=1
)
ax6.set_xlim(years.min(), years.max())
ax6.set_ylim(-3, 3)
ax6.set_ylabel('Standardized', fontsize=12)
ax6.set_xlabel('Year (DJF Season)', fontsize=12)
ax6.set_title('EOF 3 4.09%', fontsize=14)
ax6.grid(True, linestyle=':', alpha=0.7)
ax6.legend(loc='upper right')


# --- 4. Final Adjustments and Saving ---
# Use tight_layout to prevent labels from overlapping
fig.tight_layout(pad=3.0) 

# Save the entire figure to a single file
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/eof/DJF/IO_DJF_19592022_ALL_MODES.png', dpi=500, bbox_inches='tight')
plt.show()