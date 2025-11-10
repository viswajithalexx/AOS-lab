#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 21:50:58 2025

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

season_name = 'ON'
#%% region specifiying
socat = pco2_1.sel(lat = slice(-30,30),lon = slice(30,110))
#%% seasons selecting
da = socat

da_on = da.sel(mtime=da.mtime.dt.month.isin([10,11]))
da_on_yearly_mean = da_on.resample(mtime='Y').mean(dim='mtime')
on_climatology = da_on_yearly_mean.mean(dim='mtime')
on_ano = da_on_yearly_mean - on_climatology
#%% detrending
import scipy.signal as sig
import time

on_ano = on_ano.fillna(0)
on_detrend = sig.detrend(on_ano, axis=0)
on_detrend[on_detrend == 0] = np.nan
on_detrend = on_detrend[2:-2]
#%% weightage
import eofs
from eofs.standard import Eof

n_time, n_lat, n_lon = on_detrend.shape
eof_weights_lat = np.sqrt(np.cos(np.deg2rad(on_ano['lat'].values)))
weights_2d_grid = np.broadcast_to(eof_weights_lat[:, np.newaxis], (n_lat, n_lon)) # takes the latitude weight and distribute it lat*lon                                                                               # all the weights now is in 30*40 [w30 rows, 40 columns]
weights_1d_space = weights_2d_grid.flatten()

on_2d = on_detrend.reshape(n_time, n_lat * n_lon)
#%% solver
solver_on = Eof(on_2d, weights=weights_1d_space)
#%% eof maps
neofs = 3
covmaps = solver_on.eofsAsCovariance(neofs=neofs)
#%% eof maps diff. modes
covmaps_mode1 = covmaps[0].reshape(n_lat,n_lon)
covmaps_mode2 = covmaps[1].reshape(n_lat,n_lon)
covmaps_mode3 = covmaps[2].reshape(n_lat,n_lon)
#%%
a,b = np.nanpercentile(covmaps_mode1,[5,95])
#%% variance
eofvar = solver_on.varianceFraction(neigs=neofs) * 100
eofvar
#%% pcs
neofs = 3
pcs = solver_on.pcs(pcscaling=1, npcs=neofs)
#%% eof1
plt.figure(figsize =(15,10),dpi = 300)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im = plt.pcolormesh(da.lon, da.lat, covmaps_mode1,vmin = -8,vmax =8,cmap=cmaps.cmp_b2r)
ax.add_feature(cf.LAND, color='grey', zorder=12)
ax.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
cbar = plt.colorbar(im, orientation='vertical', shrink=0.7)
cbar.set_label('$\mu$atm',labelpad = 5)

# Setting the suggested title label
ax.set_title(f'IO:{season_name}:1959-2022', fontsize=16,pad =20)
# ➡️ Add the Subtitle for Variance Explained using ax.text() ⬅️
ax.text(
    1.09, 
    0.92, 
    f'{eofvar[0]:.2f} %',  # <-- CORRECTED: Formatted into one string
    horizontalalignment='center',
    transform=ax.transAxes,
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
plt.savefig(f'/home/bobco-08/24cl05012/CO2/plot/eof/{season_name}/IO_{season_name}_19592022_eof1.png', dpi=500, bbox_inches='tight')
plt.show()

#%%  PCs1

years = np.arange(1959,2023,1)

plt.figure(figsize=(15,6), dpi=500) 

plt.plot(years, pcs[:, 0], linewidth=2.5, color='darkblue', label='PC1')

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
plt.ylim(-4,4)
plt.ylabel('Standardized', fontsize=12)
plt.xlabel('Year', fontsize=12)
plt.title(f'EOF 1:{season_name}:{eofvar[0]:.2f}%', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig(f'/home/bobco-08/24cl05012/CO2/plot/eof/{season_name}/IO_{season_name}_19592022_pcs1.pdf', dpi=500, bbox_inches='tight')
plt.show()

#%% eof2
plt.figure(figsize =(15,10),dpi = 300)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im = plt.pcolormesh(da.lon, da.lat, covmaps_mode2,vmin = -8,vmax =8,cmap=cmaps.cmp_b2r)
ax.add_feature(cf.LAND, color='grey', zorder=12)
ax.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
cbar = plt.colorbar(im, orientation='vertical', shrink=0.7)
cbar.set_label('$\mu$atm',labelpad = 5)

# Setting the suggested title label
ax.set_title(f'IO:{season_name}:1959-2022', fontsize=16,pad =20)
# ➡️ Add the Subtitle for Variance Explained using ax.text() ⬅️
ax.text(
    1.09,                                         # Center position (0.5 is the center of the x-axis)
    0.92,                                        # Vertical position (0.95 is just below the top edge)
    f'{eofvar[1]:.2f} %', # The text string
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
plt.savefig(f'/home/bobco-08/24cl05012/CO2/plot/eof/{season_name}/IO_{season_name}_19592022_eof2.png', dpi=500, bbox_inches='tight')
plt.show()
#%% pcs2
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
plt.ylim(-4,4)
plt.ylabel('Standardized', fontsize=12)
plt.xlabel('Year', fontsize=12)
plt.title(f'EOF 2:{season_name}:{eofvar[1]:.2f}%', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig(f'/home/bobco-08/24cl05012/CO2/plot/eof/{season_name}/IO_{season_name}_19592022_pcs2.pdf', dpi=500, bbox_inches='tight')
plt.show()
#%% eof3
plt.figure(figsize =(15,10),dpi = 300)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im = plt.pcolormesh(da.lon, da.lat, covmaps_mode3,vmin = -8,vmax =8,cmap=cmaps.cmp_b2r)
ax.add_feature(cf.LAND, color='grey', zorder=12)
ax.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
cbar = plt.colorbar(im, orientation='vertical', shrink=0.7)
cbar.set_label('$\mu$atm',labelpad = 5)

# Setting the suggested title label
ax.set_title(f'IO:{season_name}:1959-2022', fontsize=16,pad =20)
# ➡️ Add the Subtitle for Variance Explained using ax.text() ⬅️
ax.text(
    1.09,                                         # Center position (0.5 is the center of the x-axis)
    0.92,                                        # Vertical position (0.95 is just below the top edge)
    f'{eofvar[2]:.2f} %', # The text string
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
plt.savefig(f'/home/bobco-08/24cl05012/CO2/plot/eof/{season_name}/IO_{season_name}_19592022_eof3.png', dpi=500, bbox_inches='tight')
plt.show()
#%% pcs3
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
plt.ylim(-4,5)
plt.ylabel('Standardized', fontsize=12)
plt.xlabel('Year', fontsize=12)
plt.title(f'EOF 3:{season_name}:{eofvar[2]:.2f}%', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig(f'/home/bobco-08/24cl05012/CO2/plot/eof/{season_name}/IO_{season_name}_19592022_pcs3.pdf', dpi=500, bbox_inches='tight')
plt.show()
#%% combined
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmaps


season_name = 'ON'
years = np.arange(1959, 2023, 1)


fig = plt.figure(figsize=(30, 30), dpi=300)

# ---------------------------------
# --- Row 1: EOF 1 and PC 1 ---
# ---------------------------------

# --- AX 1: EOF 1 Map (Top-Left) ---
ax1 = fig.add_subplot(3, 2, 1, projection=ccrs.PlateCarree())
ax1.coastlines(zorder=15)
im1 = ax1.pcolormesh(da.lon, da.lat, covmaps_mode1, vmin=-8, vmax=8, cmap=cmaps.cmp_b2r, transform=ccrs.PlateCarree())
ax1.add_feature(cf.LAND, color='grey', zorder=12)
gl1 = ax1.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl1.top_labels = False
gl1.right_labels = False
gl1.xlabel_style = {'fontsize': 26}
gl1.ylabel_style = {'fontsize': 26}
cbar1 = plt.colorbar(im1, ax=ax1, orientation='vertical', shrink=0.7)
cbar1.ax.tick_params(labelsize=26)
cbar1.set_label('$\mu$atm', labelpad=5, fontsize=28)
ax1.set_title(f'IO:{season_name}:1959-2022', fontsize=28, pad=20)
ax1.text(
    1.09, 0.90, f'{eofvar[0]:.2f} %',
    horizontalalignment='center',
    transform=ax1.transAxes,
    fontsize=28,
    color='k'
)
ax1.text(
    0.024, 1.02, 'EOF 1',
    horizontalalignment='center',
    transform=ax1.transAxes,
    fontsize=28,
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
ax2.set_ylim(-4, 4)
ax2.set_ylabel('Standardized', fontsize=28)
ax2.set_xlabel('Year', fontsize=28)
ax2.tick_params(axis='both', which='major', labelsize=26)
ax2.set_title(f'EOF 1 ({eofvar[0]:.2f} %)', fontsize=28)
ax2.grid(True, linestyle=':', alpha=0.7)
ax2.legend(loc='upper right', fontsize=28)

# ADD SOLID BOUNDARIES TO AX2
for spine in ax2.spines.values():
    spine.set_visible(True)
    spine.set_edgecolor('black')
    spine.set_linewidth(2)

# ---------------------------------
# --- Row 2: EOF 2 and PC 2 ---
# ---------------------------------

# --- AX 3: EOF 2 Map (Middle-Left) ---
ax3 = fig.add_subplot(3, 2, 3, projection=ccrs.PlateCarree())
ax3.coastlines(zorder=15)
im3 = ax3.pcolormesh(da.lon, da.lat, covmaps_mode2, vmin=-8, vmax=8, cmap=cmaps.cmp_b2r, transform=ccrs.PlateCarree())
ax3.add_feature(cf.LAND, color='grey', zorder=12)
gl3 = ax3.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl3.top_labels = False
gl3.right_labels = False
gl3.xlabel_style = {'fontsize': 26}
gl3.ylabel_style = {'fontsize': 26}
cbar3 = plt.colorbar(im3, ax=ax3, orientation='vertical', shrink=0.7)
cbar3.ax.tick_params(labelsize=26)
cbar3.set_label('$\mu$atm', labelpad=5, fontsize=28)
ax3.set_title(f'IO:{season_name}:1959-2022', fontsize=28, pad=20)
ax3.text(
    1.09, 0.90, f'{eofvar[1]:.2f} %',
    horizontalalignment='center',
    transform=ax3.transAxes,
    fontsize=28,
    color='k'
)
ax3.text(
    0.024, 1.02, 'EOF 2',
    horizontalalignment='center',
    transform=ax3.transAxes,
    fontsize=28,
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
ax4.set_ylim(-4, 4)
ax4.set_ylabel('Standardized', fontsize=28)
ax4.set_xlabel('Year', fontsize=28)
ax4.tick_params(axis='both', which='major', labelsize=26)
ax4.set_title(f'EOF 2 ({eofvar[1]:.2f} %)', fontsize=28)
ax4.grid(True, linestyle=':', alpha=0.7)
ax4.legend(loc='upper right', fontsize=28)

# ADD SOLID BOUNDARIES TO AX4
for spine in ax4.spines.values():
    spine.set_visible(True)
    spine.set_edgecolor('black')
    spine.set_linewidth(2)

# ---------------------------------
# --- Row 3: EOF 3 and PC 3 ---
# ---------------------------------

# --- AX 5: EOF 3 Map (Bottom-Left) ---
ax5 = fig.add_subplot(3, 2, 5, projection=ccrs.PlateCarree())
ax5.coastlines(zorder=15)
im5 = ax5.pcolormesh(da.lon, da.lat, covmaps_mode3, vmin=-8, vmax=8, cmap=cmaps.cmp_b2r, transform=ccrs.PlateCarree())
ax5.add_feature(cf.LAND, color='grey', zorder=12)
gl5 = ax5.gridlines(visible=False, draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl5.top_labels = False
gl5.right_labels = False
gl5.xlabel_style = {'fontsize': 26}
gl5.ylabel_style = {'fontsize': 26}
cbar5 = plt.colorbar(im5, ax=ax5, orientation='vertical', shrink=0.7)
cbar5.ax.tick_params(labelsize=26)
cbar5.set_label('$\mu$atm', labelpad=5, fontsize=28)
ax5.set_title(f'IO:{season_name}:1959-2022', fontsize=28, pad=20)
ax5.text(
    1.09, 0.90, f'{eofvar[2]:.2f} %',
    horizontalalignment='center',
    transform=ax5.transAxes,
    fontsize=28,
    color='k'
)
ax5.text(
    0.024, 1.02, 'EOF 3',
    horizontalalignment='center',
    transform=ax5.transAxes,
    fontsize=28,
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
ax6.set_ylim(-4, 5)
ax6.set_ylabel('Standardized', fontsize=28)
ax6.set_xlabel('Year', fontsize=28)
ax6.tick_params(axis='both', which='major', labelsize=26)
ax6.set_title(f'EOF 3 ({eofvar[2]:.2f} %)', fontsize=28)
ax6.grid(True, linestyle=':', alpha=0.7)
ax6.legend(loc='upper right', fontsize=28)

# ADD SOLID BOUNDARIES TO AX6
for spine in ax6.spines.values():
    spine.set_visible(True)
    spine.set_edgecolor('black')
    spine.set_linewidth(2)


# --- 4. Final Adjustments and Saving ---
fig.tight_layout(pad=3.0)

# Save the entire figure to a single file

plt.savefig(f'/home/bobco-08/24cl05012/CO2/plot/eof/{season_name}/IO_{season_name}_19592022_ALL_MODES.png', dpi=500, bbox_inches='tight')
plt.show()
#%%
exact_yr = covmaps_mode1*pcs[41,1]

plt.figure()
plt.figure(figsize =(15,10),dpi = 300)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im = plt.pcolormesh(da.lon, da.lat,exact_yr, vmin = -8,vmax =8,cmap=cmaps.cmp_b2r)
ax.add_feature(cf.LAND, color='grey', zorder=12)
ax.gridlines(visible=True, draw_labels=True, dms=True, x_inline=False, y_inline=False)
cbar = plt.colorbar(im, orientation='vertical', shrink=0.7)
cbar.set_label('EOF-1 Spatial Pattern ($\mu$atm)',labelpad =20)