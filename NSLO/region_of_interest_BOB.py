# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 19:38:35 2025

@author: user
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Define the region of interest
lon_min, lon_max = 82, 90
lat_min, lat_max = 8, 13

# Create the figure and axis
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8, 8))
ax.set_extent([75, 95, 5, 25], crs=ccrs.PlateCarree())

# Add features like coastlines and land
ax.add_feature(cfeature.LAND, edgecolor='black', zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=0.8, zorder=2)

# Draw the rectangle for the region of interest
ax.plot(
    [lon_min, lon_min, lon_max, lon_max, lon_min],
    [lat_min, lat_max, lat_max, lat_min, lat_min],
    color='red', linewidth=2, transform=ccrs.PlateCarree(), zorder=3
)

# Add gridlines
gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linewidth=0.5, color='gray', linestyle='--')
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()

# Title and show
ax.set_title('Region of Interest over the Bay of Bengal (BOB) (13N-8N,82E-90E)', fontsize=12)
plt.show() 