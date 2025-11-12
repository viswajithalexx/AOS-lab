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



file = "/home/bobco-08/24cl05012/CO2/data/oc_v2025.pCO2.nc"

data = xr.open_dataset(file,decode_timedelta=True)


lat =  data['lat']
lon = data['lon']
time = data['mtime']
co2_ocean = data['pCO2']

co2_flux = co2_ocean.sel(lat = slice(-30,30),lon = slice(30,110))


co2mon = co2_flux.groupby('mtime.month').mean()
co2mon = co2mon.where(co2mon != 0)


#%%

fig, axs = plt.subplots(3, 4, figsize=(18, 14), subplot_kw={'projection': ccrs.PlateCarree()})

months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

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
    gl.xlabel_style = {'fontsize': 22}
    gl.ylabel_style = {'fontsize': 22}
cbar_ax = fig.add_axes([1.02, 0.30, 0.015, 0.6]) 

# 2. Draw the colorbar in that new axis
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label(label='$\mu$atm', labelpad=10, fontsize = 26)
# 4. Corrected main title to match your data
plt.suptitle('Monthly Mean pCO₂ in Indian Ocean (1957–2024)', fontsize=16)

# 5. FIXED: Changed bottom rect from 0.3 to 0.03 to remove the gap
plt.tight_layout(rect=[0, 0.2, 1.0, 1.01])
plt.savefig('/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/monthly_mean_pco2_IO.png', dpi=500, bbox_inches='tight') 
plt.show()
#%%

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmaps
import imageio.v2 as imageio
import os

# --- Configuration ---
months_labels = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
output_dir = '/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/monthly_animation'
gif_filename = os.path.join(output_dir, 'monthly_mean_pco2_IO_animation.gif')

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# ----------------------------------------------------
# --- DUMMY DATA: REPLACE WITH YOUR ACTUAL DATA LOAD ---
# ----------------------------------------------------
# Replace these lines with your actual data loading:
# e.g., co2_flux = xr.open_dataset('your_data_file.nc')
# e.g., co2mon = co2_flux['pco2'].groupby('time.month').mean().values
lon_dummy = np.linspace(40, 100, 72)
lat_dummy = np.linspace(-30, 30, 36)
co2_flux = type('obj', (object,), {'lon': lon_dummy, 'lat': lat_dummy})()
co2mon = [np.random.rand(len(lat_dummy), len(lon_dummy)) * 180 + 300 for _ in range(12)]
# ----------------------------------------------------
# ----------------------------------------------------

# Plot levels and configuration
plot_levels = np.arange(300, 480, 10)
temp_image_files = [] # List to store paths of temporary image files

print(f"Generating 12 monthly frames for GIF in: {output_dir}")

for i in range(12): # Loop through each month
    fig = plt.figure(figsize=(18, 14), dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # Plot the monthly mean pCO2 data
    im = ax.contourf(co2_flux.lon, co2_flux.lat, co2mon[i],
                     levels=plot_levels,
                     cmap=cmaps.cmp_b2r,
                     extend='both',
                     transform=ccrs.PlateCarree())

    # Title specific to the month
    ax.set_title(f'Monthly Mean pCO$_2$ in Indian Ocean ({months_labels[i]})', fontsize=28, pad=20)

    # Map features
    ax.coastlines(zorder=15)
    ax.add_feature(cfeature.LAND, color='grey', zorder=5)
    ax.set_extent([40, 100, -30, 30], crs=ccrs.PlateCarree()) # Indian Ocean extent

    # Gridlines and labels (set to 26)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontsize': 26}
    gl.ylabel_style = {'fontsize': 26}

    # Add solid black border (linewidth 2) to the map box
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor('black')
        spine.set_linewidth(2)

    # Colorbar (label and ticks set to 28/26)
    cbar = fig.colorbar(im, ax=ax, orientation='vertical', shrink=0.7)
    cbar.ax.tick_params(labelsize=26)
    cbar.set_label('pCO$_2$ ($\mu$atm)', labelpad=10, fontsize=28)

    plt.tight_layout(rect=[0, 0.03, 1.0, 0.95])

    # Save frame and clean up
    temp_filename = os.path.join(output_dir, f'frame_{i:02d}.png')
    plt.savefig(temp_filename, dpi=150, bbox_inches='tight')
    temp_image_files.append(temp_filename)
    plt.close(fig)

# --- Create the GIF ---
print("\nCreating GIF...")
with imageio.get_writer(gif_filename, mode='I', fps=1) as writer: # 1 frame per second (fps=1)
    for filename in temp_image_files:
        image = imageio.imread(filename)
        writer.append_data(image)
        os.remove(filename) # Clean up temporary image files

print(f"GIF successfully saved to: {gif_filename}")
#%%
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmaps
import os # Import os for file path management

# --- Configuration ---
# Full month names for titles and file saving
full_months = ['January', 'February', 'March', 'April', 'May', 'June',
               'July', 'August', 'September', 'October', 'November', 'December']
# Abbreviated month names for file saving (if preferred)
abbr_months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# 1. Your plot levels
plot_levels = np.arange(300, 480, 10)

# Define the base output path
base_output_dir = '/home/bobco-08/24cl05012/CO2/plot/sem_3 plots/monthly_individual_pco2'
os.makedirs(base_output_dir, exist_ok=True) # Ensure the output directory exists

# ----------------------------------------------------
# --- DUMMY DATA PLACEHOLDERS ---
# NOTE: You MUST have 'co2_flux' (with .lon, .lat) and 'co2mon' (a list/array of 12 maps) defined
# in your actual environment for this code to run.
# For example, if you ran the previous animation code, co2_flux and co2mon are needed.
# ----------------------------------------------------

for i in range(12):
    # Create a new figure and a single subplot for each month
    fig = plt.figure(figsize=(9, 7), dpi=500) # Smaller, single-panel size
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # --- Plotting the Data ---
    im = ax.contourf(co2_flux.lon, co2_flux.lat, co2mon[i],
                     levels=plot_levels,
                     cmap=cmaps.cmp_b2r,
                     extend='both',
                     transform=ccrs.PlateCarree())

    # --- Title and Month Label ---
    month_name = full_months[i]
    abbr_name = abbr_months[i]

    ax.set_title(f'{month_name} Mean pCO$_2$ in Indian Ocean (1957–2024)', fontsize=16, pad=10)

    # Place the abbreviated month name on the plot (optional, but good for quick ID)
    ax.text(0.1, 0.95, abbr_name,
            transform=ax.transAxes,
            fontsize=16,
            fontweight='bold',
            color='black',
            zorder=20,
            ha='center'
            )

    # --- Map Features and Gridlines ---
    ax.coastlines(zorder=15)
    ax.add_feature(cfeature.LAND, color='grey', zorder=5)
    # Adjust extent to Indian Ocean region (adjust as needed)
    ax.set_extent([40, 100, -30, 30], crs=ccrs.PlateCarree())

    # Gridlines and Labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontsize': 10}
    gl.ylabel_style = {'fontsize': 10}
    
    # --- Colorbar ---
    cbar = fig.colorbar(im, ax=ax, orientation='vertical', shrink=0.8)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('pCO$_2$ ($\mu$atm)', labelpad=5, fontsize=12)

    # --- Final Layout and Saving ---
    plt.tight_layout(rect=[0, 0.0, 1.0, 0.95])
    
    # Save the individual plot using the abbreviated month name in the filename
    output_filename = os.path.join(base_output_dir, f'monthly_mean_pco2_IO_{abbr_name}.tiff')
    plt.savefig(output_filename, dpi=500, bbox_inches='tight')
    plt.close(fig) # Close the figure to free up memory

print(f"Successfully generated 12 individual plots in: {base_output_dir}")