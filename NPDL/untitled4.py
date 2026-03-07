#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 19:35:53 2026

@author: bobco-08
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Rectangle

# -----------------------------
# Define regions as dictionaries
# -----------------------------
regions = [
    {"name": "Indian Ocean", "lat": (-30, 30), "lon": (35, 110), "color": "red"},
    {"name": "Labrador Sea", "lat": (53, 63), "lon": (302.08, 317.2), "color": "blue"},
    {"name": "Antarctic", "lat": (-65, -59), "lon": (261.73, 288.28), "color": "green"}
]

# -----------------------------
# Create global map
# -----------------------------
fig = plt.figure(figsize=(14, 7), dpi=300)
ax = plt.axes(projection=ccrs.PlateCarree())

# Add map features
ax.set_global()
ax.coastlines(linewidth=0.8)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

# -----------------------------
# Plot boxes for each region
# -----------------------------
for reg in regions:
    lon_min, lon_max = reg["lon"]
    lat_min, lat_max = reg["lat"]

    # Draw rectangle
    box = Rectangle(
        (lon_min, lat_min),             # lower-left corner
        lon_max - lon_min,              # width
        lat_max - lat_min,              # height
        linewidth=2,
        edgecolor=reg["color"],
        facecolor='none',
        zorder=5,
        transform=ccrs.PlateCarree()
    )
    ax.add_patch(box)
    
    # Add label at top-left corner
    ax.text(
        lon_min, lat_max + 1,
        reg["name"],
        color=reg["color"],
        fontsize=10,
        weight='bold',
        transform=ccrs.PlateCarree()
    )

plt.title("Global Regions of Interest", fontsize=14, weight='bold')
plt.show()
