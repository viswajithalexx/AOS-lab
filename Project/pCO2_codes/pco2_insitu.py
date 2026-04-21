#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:54:20 2026

@author: bobco-08
"""

import pandas as pd

# --------------------------------------------------
# 1. Read the original SOCAT Excel file
# --------------------------------------------------
file_path = "/home/bobco-08/Desktop/socat_data_edited.xlsx"   # change to your file path
df = pd.read_excel(file_path)

# --------------------------------------------------
# 2. Define Expocodes for each region
# --------------------------------------------------
#%%
EIO_expocodes = [
"06BE20140723","06MT19950714","316N19950124","316N19950310",
"316N19950715","316N19950829","316N19951205","325020240408",
"33MW19950321","33MW19950427","33MW19950531","33MW19950712",
"33MW19950731","33MW19950902","33MW19951028","33RO19990222",
"33RO19990304","33RO19990326","33RO19990430","33RO20180423",
"33RO20180519","33RO20180616","33RR20070322","33RR20160321",
"35MF19991018","49NZ20010922","76XL20080304"
]

NWIO_expocodes = [
"06BE19970515","06BE19970612","06MT19950714","316N19950715",
"316N19950829","325019950109","325019950210","325019950315",
"325019950721","325019951101","325019951204","33MW19950321",
"33MW19950427","33MW19950531","33MW19950712","33MW19950731",
"33RO20180519","33RO20180616","35MF19991018"
]

NAS_expocodes = [
"06BE19970515","06BE19970612","06MT19950714","316N19950715",
"325019950109","325019950721","325019951101","33MW19950427",
"33MW19950531","33MW19950712","33MW19950731"
]

ESIO_expocodes = [
"316N19950124","316N19950310","316N19950829","316N19951205",
"325020240408","33MW19950902","33MW19951028","33RO19990430",
"33RO20180616","33RR20070322","33RR20160321","49EM19950123",
"49HM19960120","49NZ20010922","76XL20080304"
]
#%%
# --------------------------------------------------
# 3. Extract rows belonging to those Expocodes
# --------------------------------------------------

EIO_data = df[df["Expocode"].isin(EIO_expocodes)]
NWIO_data = df[df["Expocode"].isin(NWIO_expocodes)]
NAS_data = df[df["Expocode"].isin(NAS_expocodes)]
ESIO_data = df[df["Expocode"].isin(ESIO_expocodes)]
#%%
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def plot_tracks_with_box(data, region_name, lon_min, lon_max, lat_min, lat_max):

    fig = plt.figure(figsize=(10,8))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Wider extent so full tracks appear
    ax.set_extent([30,120,-30,30], crs=ccrs.PlateCarree())

    # Land and coastlines
    ax.add_feature(cfeature.LAND, facecolor='lightgrey')
    ax.add_feature(cfeature.COASTLINE, linewidth=1)

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False

    # Plot cruise tracks
    for expocode, group in data.groupby("Expocode"):
        ax.scatter(group["longitude [dec.deg.E]"],
               group["latitude [dec.deg.N]"],
               s=6,
               transform=ccrs.PlateCarree(),
               label=expocode)

    # Draw region box
    box_lons = [lon_min, lon_max, lon_max, lon_min, lon_min]
    box_lats = [lat_min, lat_min, lat_max, lat_max, lat_min]

    ax.plot(box_lons, box_lats,
            color='red',
            linewidth=2,
            transform=ccrs.PlateCarree())

    plt.title(f"Cruise Tracks with Region Box ({region_name})", fontsize=15)

    plt.legend(fontsize=7, bbox_to_anchor=(1.05,1), loc='upper left')

    plt.tight_layout()

    plt.show()

#%%
# --------------------------------------------------
# Plot regions
# --------------------------------------------------

plot_tracks_with_box(EIO_data,
                     "EIO",
                     49,92,
                     -6.5,5)

plot_tracks_with_box(NWIO_data,
                     "NWIO",
                     45,65,
                     5,22.5)

plot_tracks_with_box(NAS_data,
                     "NAS",
                     56,70,
                     22.5,28)

plot_tracks_with_box(ESIO_data,
                     "ESIO",
                     92,109,
                     -6.6,8)

#%%


def plot_intersecting_tracks_scatter(df, region_name, lon_min, lon_max, lat_min, lat_max):

    fig = plt.figure(figsize=(10,8))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Large map extent so full cruise path appears
    ax.set_extent([30,120,-30,30], crs=ccrs.PlateCarree())

    # Land and coastlines
    ax.add_feature(cfeature.LAND, facecolor='lightgrey')
    ax.add_feature(cfeature.COASTLINE)

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False

    intersecting_expocodes = []

    for expocode, group in df.groupby("Expocode"):

        # Check if cruise intersects region
        inside = group[
            (group["longitude [dec.deg.E]"] >= lon_min) &
            (group["longitude [dec.deg.E]"] <= lon_max) &
            (group["latitude [dec.deg.N]"] >= lat_min) &
            (group["latitude [dec.deg.N]"] <= lat_max)
        ]

        if len(inside) > 0:

            intersecting_expocodes.append(expocode)

            # Scatter plot of ALL observations of that cruise
            ax.scatter(group["longitude [dec.deg.E]"],
                       group["latitude [dec.deg.N]"],
                       s=6,
                       transform=ccrs.PlateCarree(),
                       label=expocode)

    # Draw region box
    box_lon = [lon_min, lon_max, lon_max, lon_min, lon_min]
    box_lat = [lat_min, lat_min, lat_max, lat_max, lat_min]

    ax.plot(box_lon, box_lat, color='red', linewidth=2, transform=ccrs.PlateCarree())

    plt.title(f"Cruise Observations Intersecting {region_name}")

    plt.legend(fontsize=7, bbox_to_anchor=(1.05,1), loc='upper left')

    plt.tight_layout()
    plt.show()

    return intersecting_expocodes

#%%

NWIO_tracks = plot_intersecting_tracks_scatter(
    df,
    "NWIO",
    45,65,
    5,22.5
)



EIO_tracks = plot_intersecting_tracks_scatter(
    df,
    "EIO",
    49,92,
    -6.5,5
)



NAS_tracks = plot_intersecting_tracks_scatter(
    df,
    "NAS",
    56,70,
    22.5,24
)

ESIO_tracks = plot_intersecting_tracks_scatter(
    df,
    "ESIO",
    92,109,
    -6.6,8
)


#%%
# --------------------------------------------------
# 4. Save extracted datasets
# --------------------------------------------------

EIO_data.to_excel("EIO_insitu_data.xlsx", index=False)
NWIO_data.to_excel("NWIO_insitu_data.xlsx", index=False)
NAS_data.to_excel("NAS_insitu_data.xlsx", index=False)
ESIO_data.to_excel("ESIO_insitu_data.xlsx", index=False)

print("Extraction complete!")

#%%
import pandas as pd

# --------------------------------------------------
# 1. Read SOCAT Excel file
# --------------------------------------------------

file_path = "/home/bobco-08/Desktop/socat_data_edited.xlsx"
df = pd.read_excel(file_path)

# Short variable names
lon = df["longitude [dec.deg.E]"]
lat = df["latitude [dec.deg.N]"]

# --------------------------------------------------
# 2. Define Expocodes for each region
# --------------------------------------------------

EIO_expocodes = [
"06BE20140723","06MT19950714","316N19950124","316N19950310",
"316N19950715","316N19950829","316N19951205","325020240408",
"33MW19950321","33MW19950427","33MW19950531","33MW19950712",
"33MW19950731","33MW19950902","33MW19951028","33RO19990222",
"33RO19990304","33RO19990326","33RO19990430","33RO20180423",
"33RO20180519","33RO20180616","33RR20070322","33RR20160321",
"35MF19991018","49NZ20010922","76XL20080304"
]

NWIO_expocodes = [
"06BE19970515","06BE19970612","06MT19950714","316N19950715",
"316N19950829","325019950109","325019950210","325019950315",
"325019950721","325019951101","325019951204","33MW19950321",
"33MW19950427","33MW19950531","33MW19950712","33MW19950731",
"33RO20180519","33RO20180616","35MF19991018"
]

NAS_expocodes = [
"06BE19970515","06BE19970612","06MT19950714","316N19950715",
"325019950109","325019950721","325019951101","33MW19950427",
"33MW19950531","33MW19950712","33MW19950731"
]

ESIO_expocodes = [
"316N19950124","316N19950310","316N19950829","316N19951205",
"325020240408","33MW19950902","33MW19951028","33RO19990430",
"33RO20180616","33RR20070322","33RR20160321","49EM19950123",
"49HM19960120","49NZ20010922","76XL20080304"
]

# --------------------------------------------------
# 3. Define regional boxes
# --------------------------------------------------

EIO_box = (
    (lon >= 49) & (lon <= 92) &
    (lat >= -6.5) & (lat <= 5)
)

NWIO_box = (
    (lon >= 45) & (lon <= 65) &
    (lat >= 5) & (lat <= 22.5)
)

NAS_box = (
    (lon >= 56) & (lon <= 70) &
    (lat >= 22.5) & (lat <= 28)
)

ESIO_box = (
    (lon >= 92) & (lon <= 109) &
    (lat >= -6.6) & (lat <= 8)
)

# --------------------------------------------------
# 4. Filter by Expocode AND region
# --------------------------------------------------

EIO_data = df[
    (df["Expocode"].isin(EIO_expocodes)) &
    EIO_box
]

NWIO_data = df[
    (df["Expocode"].isin(NWIO_expocodes)) &
    NWIO_box
]

NAS_data = df[
    (df["Expocode"].isin(NAS_expocodes)) &
    NAS_box
]

ESIO_data = df[
    (df["Expocode"].isin(ESIO_expocodes)) &
    ESIO_box
]

# --------------------------------------------------
# 5. Print number of observations
# --------------------------------------------------

print("EIO points:", len(EIO_data))
print("NWIO points:", len(NWIO_data))
print("NAS points:", len(NAS_data))
print("ESIO points:", len(ESIO_data))

# --------------------------------------------------
# 6. Save filtered datasets
# --------------------------------------------------

EIO_data.to_csv("/home/bobco-08/Desktop/EIO_insitu_points.csv", index=False)
NWIO_data.to_csv("/home/bobco-08/Desktop/NWIO_insitu_points.csv", index=False)
NAS_data.to_csv("/home/bobco-08/Desktop/NAS_insitu_points.csv", index=False)
ESIO_data.to_csv("/home/bobco-08/Desktop/ESIO_insitu_points.csv", index=False)

print("Files saved successfully.")