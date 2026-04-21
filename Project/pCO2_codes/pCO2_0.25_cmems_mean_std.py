#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:53:28 2026

@author: bobco-08
"""



# importing libraries

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmaps


file = "/home/bobco-08/24cl05012/CO2/data/data_1/cmems/bio_var/cmems_obs-mob_glo_bgc-car_my_irr-i_1775027250539_30E-119.88E_29.88S-29.88N.nc"

data = xr.open_dataset(file)

lat =  data['latitude']
lon = data['longitude']
time = data['time']
pco2 = data['spco2']
# surface douwnward flux of total carbon
pco2v = pco2.values
# for IO

pco2_io = pco2.sel(time = slice('1994-01-01','2024-12-01'))

#%% annual mean of spco2

# adjusting the co2f data
co2f_clim = pco2_io.mean(dim=('time'))

# plotting 
flux_levels = np.arange(300,460,5)
ftl = np.arange(300,460,10)
fig = plt.figure(figsize=(18,14))
ax = plt.axes(projection = ccrs.PlateCarree())
ax.set_extent([30, 120, -30, 30], crs=ccrs.PlateCarree())
im = plt.contourf(lon,lat,co2f_clim,cmap = cmaps.WhiteBlueGreenYellowRed,levels = flux_levels,
                  transform = ccrs.PlateCarree(),extend = 'both')
ax.coastlines(zorder = 12)
land = cf.NaturalEarthFeature(
    'physical', 'land', '10m',
    edgecolor='black',
    facecolor='gray'
)
ax.add_feature(land, zorder=2)


# contour
cs = plt.contour(lon, lat, co2f_clim,levels=flux_levels,colors='k',
                 lineslinewidths=0.5,transform=ccrs.PlateCarree())
plt.clabel(cs, fmt='%d', fontsize=17)


# ticks and gridlines
xticks = np.arange(35,120,10)
yticks = np.arange(-30,40,10)
ax.set_xticks(xticks, crs=ccrs.PlateCarree()) # for adding ticks 
ax.set_yticks(yticks, crs=ccrs.PlateCarree())
ax.tick_params(axis='both',
               direction='out',   # makes arrow-like ticks
               length=5,
               width=1.2,
               labelsize=20)
gl = ax.gridlines(draw_labels=False, linewidth=0.5, alpha=0.1)
gl.xlocator = plt.FixedLocator(xticks)
gl.ylocator = plt.FixedLocator(yticks)

gl.top_labels = False
gl.right_labels = False
# plt.setp(ax.get_xticklabels(), fontweight='bold')
# plt.setp(ax.get_yticklabels(), fontweight='bold')

#colorbar
cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.70]) #manual colorbar
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical',ticks = ftl)
cbar.ax.tick_params(axis='both', which='major',
               direction='out',   # ticks pointing outward
               length=6,          # tick length
               width=1,labelsize=20)

#title
cbar.set_label('\u00b5atm',fontsize = 20,labelpad=12,
               fontweight = 'bold')
ax.set_title('Climatology of pCO$_2$ in Indian Ocean (1994-2024)',
             fontsize = 24,y =1.01,fontweight='bold')
plt.show()
#%% djf mean pco2

pco2_djf = pco2_io.resample(time = 'QS-DEC').mean(dim=('time'))
pco2_djf = pco2_djf.isel(time=slice(1, None))
pco2_djf = pco2_djf.isel(time = pco2_djf.time.dt.month == 12)

djf_mean = pco2_djf.mean(dim ='time')

#%% mam mean pco2

pco2_mam = pco2_io.sel(time= pco2_io.time.dt.month.isin([3,4,5]))
pco2_mam = pco2_mam.resample(time='YE').mean(dim='time')
mam_mean = pco2_mam.mean(dim='time')

#%% jjas mean pco2
pco2_jjas = pco2_io.sel(time= pco2_io.time.dt.month.isin([6,7,8,9]))
pco2_jjas = pco2_jjas.resample(time='YE').mean(dim='time')
jjas_mean = pco2_jjas.mean(dim='time')

#%% on mean pco2

pco2_on = pco2_io.sel(time= pco2_io.time.dt.month.isin([10,11]))
pco2_on = pco2_on.resample(time='YE').mean(dim='time')
on_mean = pco2_on.mean(dim='time')

#%% Seasonal climatology of pCO₂ in Indian Ocean (1980–2019)

fig, axs = plt.subplots(2, 2, figsize=(18,14),
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace':0.06,'hspace': -0.09})
season = ['DJF','MAM','JJAS',"ON"]
co2f_season = [djf_mean,mam_mean,jjas_mean,on_mean]




season_lv = np.arange(300,460,5)
slv = np.arange(300,460,10)
#plotting

for i, ax in enumerate(axs.flat):
    im = ax.contourf(lon,lat,co2f_season[i],levels =season_lv,
                     cmap= cmaps.WhiteBlueGreenYellowRed,
                     transform=ccrs.PlateCarree(),extend = 'both') 
    ax.coastlines(zorder =15)
    ax.coastlines(zorder = 12)
    ax.set_extent([30, 120, -30, 30], crs=ccrs.PlateCarree())
    land = cf.NaturalEarthFeature(
        'physical', 'land', '10m',
        edgecolor='black',
        facecolor='gray'
    )
    ax.add_feature(land, zorder=2)

# subplot individual title

    ax.text(0.75, 0.98, f"Mean ({season[i]})",transform=ax.transAxes,         
        fontsize=18, zorder= 18,fontweight='bold',
        va='top', ha='right',)       

#contours  
 
    contours = ax.contour(lon,lat,co2f_season[i],levels = 25,
                              colors='black', linewidths=0.6)
    ax.clabel(contours, inline=True, fontsize=12, fmt="%.0f")        

# ticks and gridlines

    ax.set_xticks(xticks, crs=ccrs.PlateCarree()) # for adding ticks 
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    
    ax.tick_params(axis='both',
                   direction='out',   # makes arrow-like ticks
                   length=5,
                   width=1.2,
                   labelsize=24)
    
      # Hide labels depending on subplot index
    if i in [0, 1]:          # top row
        ax.tick_params(labelbottom=False)
    
    if i in [1, 3]:          # right column
        ax.tick_params(labelleft=False)

#colorbar

cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.70]) 
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical',ticks =slv )
cbar.ax.tick_params(axis='both', which='major',
               direction='out',   # ticks pointing outward
               length=6,          # tick length
               width=1,labelsize=24)


#titles
cbar.set_label('\u00b5atm',fontsize = 24,labelpad=13,
               fontweight = 'bold')
# fig.suptitle(
#     'Seasonal Air–Sea CO$_2$ Flux in Indian Ocean (1994–2024)',
#     fontsize=24,
#     fontweight='bold',
#     y=0.88
# )
plt.show()
#%% standard deviation of seasonal plots 

djf_std = pco2_djf.std(dim =('time'))
jjas_std = pco2_jjas.std(dim =('time'))
mam_std = pco2_mam.std(dim =('time'))
on_std = pco2_on.std(dim =('time'))

fig, axs = plt.subplots(2, 2, figsize=(18,14),
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace':0.06,'hspace': -0.12})
season = ['DJF','MAM','JJAS',"ON"]
season_std = [djf_std,mam_std,jjas_std,on_std]

std_lv = np.arange(8,28,1)
stlv = np.arange(8,28,2)
#plotting

for i, ax in enumerate(axs.flat):
    im = ax.contourf(lon,lat,season_std[i],levels =std_lv,
                     cmap= cmaps.WhiteBlueGreenYellowRed,
                     transform=ccrs.PlateCarree(),extend = 'both') 
    ax.coastlines(zorder =15)
    ax.coastlines(zorder = 12)
    ax.set_extent([30, 120, -30, 30], crs=ccrs.PlateCarree())
    land = cf.NaturalEarthFeature(
        'physical', 'land', '10m',
        edgecolor='black',
        facecolor='gray'
    )
    ax.add_feature(land, zorder=2)

# subplot individual title

    ax.text(0.75, 0.98, f"Std ({season[i]})",transform=ax.transAxes,         
        fontsize=18, zorder= 18,fontweight='bold',
        va='top', ha='right',)       

# contours  
 
    contours = ax.contour(lon,lat,season_std[i],levels = 25,
                              colors='black', linewidths=0.6)
    ax.clabel(contours, inline=True, fontsize=12, fmt="%.0f")        

# ticks and gridlines

    ax.set_xticks(xticks, crs=ccrs.PlateCarree()) # for adding ticks 
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    
    ax.tick_params(axis='both',
                   direction='out',   # makes arrow-like ticks
                   length=5,
                   width=1.2,
                   labelsize=24)
    
      # Hide labels depending on subplot index
    if i in [0, 1]:          # top row
        ax.tick_params(labelbottom=False)
    
    if i in [1, 3]:          # right column
        ax.tick_params(labelleft=False)

#colorbar

cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.70]) 
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical',ticks =stlv )
cbar.ax.tick_params(axis='both', which='major',
               direction='out',   # ticks pointing outward
               length=6,          # tick length
               width=1,labelsize=24)


#titles
cbar.set_label('\u00b5atm',fontsize = 24,labelpad=12,
               fontweight = 'bold')
# fig.suptitle(
#     'Seasonal Air–Sea CO$_2$ Flux in Indian Ocean (1994–2024)',
#     fontsize=24,
#     fontweight='bold',
#     y=0.88
# )
plt.show()

#%% Standard deviation of monthly pCO₂ in Indian Ocean (1980 -2019)

mon_mean = pco2_io.groupby('time.month').mean()

mon_ano = pco2_io.groupby('time.month') - mon_mean

std_mon = mon_ano.groupby('time.month').std('time')

#%%
fig, axs = plt.subplots(3, 4, figsize=(18, 15),dpi = 400,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace': 0.08, 'hspace': -0.55})
x1 =np.arange(35,120,15)
y1 = np.arange(-30,40,15)
months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

month_lv = np.arange(300,450, 5)
mlv = np.arange(300,450,10)

for i, ax in enumerate(axs.flat):

    im = ax.contourf(lon, lat, mon_mean[i],
                     levels=month_lv,
                     cmap= cmaps.WhiteBlueGreenYellowRed,
                     extend='both',
                     transform=ccrs.PlateCarree())

    ax.coastlines(zorder=12)

    land = cf.NaturalEarthFeature(
        'physical', 'land', '10m',
        edgecolor='black',
        facecolor='gray'
    )
    ax.add_feature(land, zorder=2)

    # subplot title (cleaner alignment)
    ax.set_title(f"{months[i]}",
                 fontsize=16,
                 fontweight='bold',
                 loc='right',
                 pad=6)

    # ticks
    ax.set_xticks(x1, crs=ccrs.PlateCarree())
    ax.set_yticks(y1, crs=ccrs.PlateCarree())

    ax.tick_params(axis='both',
                   direction='out',
                   length=4,
                   width=1,
                   labelsize=14)

    # show only left & bottom ticks (clean layout)
    if i < 8:   # first 2 rows
        ax.tick_params(labelbottom=False)

    if i % 4 != 0:   # not first column
        ax.tick_params(labelleft=False)

# colorbar (better aligned)
cbar_ax = fig.add_axes([0.92, 0.21, 0.02, 0.59])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', ticks=mlv)

cbar.ax.tick_params(direction='out',
                    length=5,
                    width=1,
                    labelsize=16)

cbar.set_label('\u00b5atm',
               fontsize=18,
               labelpad=10,
               fontweight='bold')

# overall spacing
# fig.subplots_adjust(left=0.06, right=0.90, top=0.95, bottom=0.07)
# fig.suptitle(
#     'Monthly Air–Sea CO$_2$ Flux in Indian Ocean (1994–2024)',
#     fontsize=24,
#     fontweight='bold',
#     y=0.87
# )
plt.show()
#%%
fig, axs = plt.subplots(3, 4, figsize=(18, 14),dpi = 400,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace': 0.08, 'hspace': -0.55})
x1 =np.arange(35,120,15)
y1 = np.arange(-30,40,15)
months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

month_lv = np.arange(0,30,1)
mlv = np.arange(0,30,2)

for i, ax in enumerate(axs.flat):

    im = ax.contourf(lon, lat, std_mon[i],
                     levels=month_lv,
                     cmap='turbo',
                     extend = 'max',
                     transform=ccrs.PlateCarree())

    ax.coastlines(zorder=12)

    land = cf.NaturalEarthFeature(
        'physical', 'land', '10m',
        edgecolor='black',
        facecolor='gray'
    )
    ax.add_feature(land, zorder=2)

    # subplot title (cleaner alignment)
    ax.set_title(f"{months[i]}",
                 fontsize=16,
                 fontweight='bold',
                 loc='right',
                 pad=6)

    # ticks
    ax.set_xticks(x1, crs=ccrs.PlateCarree())
    ax.set_yticks(y1, crs=ccrs.PlateCarree())

    ax.tick_params(axis='both',
                   direction='out',
                   length=4,
                   width=1,
                   labelsize=14)

    # show only left & bottom ticks (clean layout)
    if i < 8:   # first 2 rows
        ax.tick_params(labelbottom=False)

    if i % 4 != 0:   # not first column
        ax.tick_params(labelleft=False)

# colorbar (better aligned)
cbar_ax = fig.add_axes([0.92, 0.21, 0.02, 0.59])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', ticks=mlv)

cbar.ax.tick_params(direction='out',
                    length=5,
                    width=1,
                    labelsize=16)

cbar.set_label('\u00b5atm',
               fontsize=18,
               labelpad=10,
               fontweight='bold')

# overall spacing
fig.subplots_adjust(left=0.06, right=0.90, top=0.95, bottom=0.07)
# fig.suptitle(
#     'Monthly Air–Sea CO$_2$ Flux in Indian Ocean (1994–2024)',
#     fontsize=24,
#     fontweight='bold',
#     y=0.87
# )
plt.show()

#%%

file ="/home/bobco-08/24cl05012/npdl/rf_monthly_1982_2020.nc"
data = xr.open_dataset(file)

data = data.rename({'TIME': 'time'})

lon = data['LONGITUDE']
lat = data['LATITUDE']
t = data['time']
rf = data['RAINFALL']

#%%