#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 07:14:08 2025

@author: bobco-08
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Define map projection
proj = ccrs.PlateCarree(central_longitude=180)

# Create figure and axis
fig = plt.figure(figsize=(10, 5))
ax = plt.axes(projection=proj)

# Set extent: [lon_min, lon_max, lat_min, lat_max]
ax.set_extent([30, 290, -30, 30], crs=ccrs.PlateCarree())

# Add coastlines and features
ax.coastlines(resolution='110m', linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

# Add gridlines and labels
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5)
gl.top_labels = False
gl.right_labels = False

# Add title
plt.title("Model Domain: Indo-Pacific Ocean (30°E–290°E, 30°S–30°N)", fontsize=12)

# Show plot
plt.show()
