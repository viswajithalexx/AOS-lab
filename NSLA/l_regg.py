# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 08:40:19 2025

@author: user
"""

import numpy as np
import xarray as xr
import matplotlib . pyplot as plt
import cartopy . feature as cf

import cartopy . crs as ccrs

# Load data
path = "C:/Users/user/Documents/24CL05012/nsla/rainfall/"
data = xr.open_mfdataset(path + "*.nc", combine='nested')

time = data ['TIME']
lon = data ['LONGITUDE']
lat = data ['LATITUDE']
rf = data ['RAINFALL']

rf = rf.sel( TIME =~(( rf . TIME . dt . month == 2) & ( rf . TIME . dt . day == 29) ) )

lon_blr = 77.35
lat_blr = 12.58
rf_blr = rf.sel( LONGITUDE = lon_blr , LATITUDE = lat_blr , method ='nearest')

rf_blr = rf_blr.fillna (0).values
fft_blr = np.fft.fft ( rf_blr )
omega = np. fft . fftfreq (len( rf_blr ) , d =1)
omega [ omega == 0] = np. nan
T_blr = 1 / omega

mask_band = ( T_blr > 20) & ( T_blr < 90)
fft_masked = np. zeros_like ( fft_blr )
fft_masked [ mask_band ] = fft_blr [ mask_band ]
rf_intra = np. fft . ifft ( fft_masked ) . real

x = xr. DataArray ( rf_intra , coords =[ rf . TIME ] , dims =[" TIME "])
x = x - x.mean ( dim =" TIME ")

rf = rf.fillna (0)

y = rf - rf.mean ( dim =" TIME ")

cov = ( x * y ).mean ( dim =" TIME ")
var_x = np. var ( x )
regression_map = cov / var_x
plt. figure ( figsize =(14 , 7) )
ax = plt. axes ( projection = ccrs . PlateCarree () )
contour = ax.contourf (
regression_map.LONGITUDE , regression_map . LATITUDE ,
regression_map , levels =21 , cmap =" coolwarm ", transform = ccrs.PlateCarree ())

ax . set_xticks (np. arange ( float ( regression_map.LONGITUDE.min () ) ,float ( regression_map.LONGITUDE.max () ) + 1 , 5) , crs = ccrs . PlateCarree () )
ax . set_yticks (np. arange ( float ( regression_map.LATITUDE.min () ) ,float ( regression_map.LATITUDE.max () ) + 1 , 5) , crs = ccrs . PlateCarree () )
ax.tick_params ( labelsize =16)
ax.set_title (" Regression : 20 90 Day Intra - Seasonal Signal vs GridRainfall ", fontsize =18)

ax.coastlines ()
ax.add_feature ( cf.BORDERS )
cbar = plt.colorbar ( contour , pad =0.05)
cbar.set_label (" Regression Coefficient ", fontsize =16)
cbar.ax.tick_params ( labelsize =14)

plt. tight_layout ()
plt. show ()