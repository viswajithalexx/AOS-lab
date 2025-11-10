#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 17:09:59 2025

@author: bobco-08
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.patches as mpatches

# --- 1. Define your regions (from your code) ---
# eio: Equatorial Indian Ocean
eio_lon_slice = slice(40, 96)
eio_lat_slice = slice(-5, 5)

# so: Somali Coast
so_lon_slice = slice(40, 70)
so_lat_slice = slice(0, 24)

# nas: Northern Arabian Sea
nas_lon_slice = slice(55, 72)
nas_lat_slice = slice(20, 27)

# ⬇️ --- NEW REGION ADDED HERE --- ⬇️
# ans: Andaman Sea
ans_lon_slice = slice(92, 110)
ans_lat_slice = slice(-5, 10)
# ⬆️ --- END OF NEW REGION --- ⬆️


# --- 2. Create the Figure and Axes ---
fig = plt.figure(figsize=(10, 8), dpi=200)
ax = plt.axes(projection=ccrs.PlateCarree())

# Set map extent to focus on the Indian Ocean
ax.set_extent([20, 120, -30, 40], crs=ccrs.PlateCarree())

# Add map features
ax.add_feature(cf.LAND, color='lightgrey', zorder=1)
ax.add_feature(cf.BORDERS, linewidth=0.5, zorder=2)
ax.coastlines(zorder=3)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False

# --- 3. Define the patches (boxes) ---
# We need: (lon_min, lat_min), width, height
# Note: We must add 'transform=ccrs.PlateCarree()' so Cartopy knows
#       these coordinates are lat/lon.

# EIO: Red, transparent
eio_box = mpatches.Rectangle(
    xy=[eio_lon_slice.start, eio_lat_slice.start],  # (lon_min, lat_min)
    width=eio_lon_slice.stop - eio_lon_slice.start,
    height=eio_lat_slice.stop - eio_lat_slice.start,
    facecolor='red',
    alpha=0.4,  # Transparency
    edgecolor='red',
    linewidth=1.5,
    transform=ccrs.PlateCarree(),
    zorder=5
)

# SO: Blue, transparent
so_box = mpatches.Rectangle(
    xy=[so_lon_slice.start, so_lat_slice.start],
    width=so_lon_slice.stop - so_lon_slice.start,
    height=so_lat_slice.stop - so_lat_slice.start,
    facecolor='blue',
    alpha=0.4,
    edgecolor='blue',
    linewidth=1.5,
    transform=ccrs.PlateCarree(),
    zorder=4  # Placed behind EIO
)

# NAS: Green, transparent
nas_box = mpatches.Rectangle(
    xy=[nas_lon_slice.start, nas_lat_slice.start],
    width=nas_lon_slice.stop - nas_lon_slice.start,
    height=nas_lat_slice.stop - nas_lat_slice.start,
    facecolor='green',
    alpha=0.4,
    edgecolor='green',
    linewidth=1.5,
    transform=ccrs.PlateCarree(),
    zorder=6
)

# ⬇️ --- NEW PATCH ADDED HERE --- ⬇️
# ANS: Purple, transparent
ans_box = mpatches.Rectangle(
    xy=[ans_lon_slice.start, ans_lat_slice.start],
    width=ans_lon_slice.stop - ans_lon_slice.start,
    height=ans_lat_slice.stop - ans_lat_slice.start,
    facecolor='purple',
    alpha=0.4,
    edgecolor='purple',
    linewidth=1.5,
    transform=ccrs.PlateCarree(),
    zorder=7
)
# ⬆️ --- END OF NEW PATCH --- ⬆️


# --- 4. Add Patches and Text to the Map ---
ax.add_patch(eio_box)
ax.add_patch(so_box)
ax.add_patch(nas_box)
ax.add_patch(ans_box) # ⬅️ Added new patch to map

# Add text labels
ax.text(eio_lon_slice.start + 2, eio_lat_slice.start + 1, 'EIO', 
        color='black', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)

ax.text(so_lon_slice.start + 2, so_lat_slice.start + 1, 'NWIO', 
        color='black', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)

ax.text(nas_lon_slice.start + 1, nas_lat_slice.start + 1, 'NAS', 
        color='black', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)

# ⬅️ Added new text to map
ax.text(ans_lon_slice.start + 1, ans_lat_slice.start + 1, 'ESIO', 
        color='black', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)

ax.set_title("pCO₂ Analysis Regions")
plt.show()