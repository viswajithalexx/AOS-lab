#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 12:14:09 2026

@author: bobco-08
"""
# importing libraries

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmaps


file = "/home/bobco-08/24cl05012/CO2/data/data_1/jena/oc_v2025.flux.nc"

data = xr.open_dataset(file)


data = data.rename({'mtime': 'time'})
data = data.rename({'lat': 'latitude'})
data = data.rename({'lon': 'longitude'})

lat =  data['latitude']
lon = data['longitude']
time = data['time']
dxyp = data['dxyp']
co2f = data['co2flux_ocean'].sel(latitude = slice(-30,30),longitude = slice(30,120),time = slice('1994-01-01','2024-12-01'))
co2f = co2f/dxyp

co2fv = co2f.values*1e15
# for IO

co2f_io = co2f*1e15
lon = co2f.longitude
lat = co2f.latitude


#%% annual mean of surface downward flux co2

# adjusting the co2f data
co2f_clim = co2f_io.mean(dim=('time'))

# plotting 
flux_levels = np.linspace(-40,40,21)
ftl = np.arange(-40,44,4)
fig = plt.figure(figsize=(18,14))
ax = plt.axes(projection = ccrs.PlateCarree())
im = plt.contourf(lon,lat,co2f_clim,cmap = cmaps.cmocean_balance,levels = flux_levels ,
                  transform = ccrs.PlateCarree(),extend = 'both')
ax.coastlines(zorder = 12)
ax.set_extent([30, 120, -29, 29], crs=ccrs.PlateCarree())
land = cf.NaturalEarthFeature(
    'physical', 'land', '10m',
    edgecolor='black',
    facecolor='gray'
)
ax.add_feature(land, zorder=2)


# contour
# cs = plt.contour(lon, lat, co2f_clim,levels=ftl,colors='k',
#                  lineslinewidths=0.5,transform=ccrs.PlateCarree())
# plt.clabel(cs, fmt='%d', fontsize=17)


# ticks and gridlines
xticks = np.arange(35,125,10)
yticks = np.arange(-20,40,10)
ax.set_xticks(xticks, crs=ccrs.PlateCarree()) # for adding ticks 
ax.set_yticks(yticks, crs=ccrs.PlateCarree())
ax.tick_params(axis='both',
               direction='out',   # makes arrow-like ticks
               length=5,
               width=1.2,
               labelsize=24)
gl = ax.gridlines(draw_labels=False, linewidth=0.5, alpha=0.1)
gl.xlocator = plt.FixedLocator(xticks)
gl.ylocator = plt.FixedLocator(yticks)

gl.top_labels = False
gl.right_labels = False
# plt.setp(ax.get_xticklabels(), fontweight='bold')
# plt.setp(ax.get_yticklabels(), fontweight='bold')

#colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.70]) #manual colorbar
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical',ticks =ftl)
cbar.ax.tick_params(axis='both', which='major',
               direction='out',   # ticks pointing outward
               length=6,          # tick length
               width=1,labelsize=24)

#title
cbar.set_label('flux (gC m$^{-2}$ yr$^{-1}$)',fontsize = 24,labelpad=12,
               fontweight = 'bold')
# ax.set_title('Air-Sea CO$_2$ flux climatology in Indian Ocean (1994-2024)\n Jena Carboscope',
#              fontsize = 24,y =1.01,fontweight='bold')
plt.show()
#%% djf mean co2_flux

co2f_djf = co2f_io.resample(time = 'QS-DEC').mean(dim=('time'))
co2f_djf = co2f_djf.isel(time=slice(1, None))
co2f_djf = co2f_djf.isel(time = co2f_djf.time.dt.month == 12)

djf_mean = co2f_djf.mean(dim ='time')

#%% mam mean co2_flux

co2f_mam = co2f_io.sel(time= co2f_io.time.dt.month.isin([3,4,5]))
co2f_mam = co2f_mam.resample(time='YE').mean(dim='time')
mam_mean = co2f_mam.mean(dim='time')

#%% jjas mean co2_flux

co2f_jjas = co2f_io.sel(time= co2f_io.time.dt.month.isin([6,7,8,9]))
co2f_jjas = co2f_jjas.resample(time='YE').mean(dim='time')
jjas_mean = co2f_jjas.mean(dim='time')

#%% on mean co2_flux

co2f_on = co2f_io.sel(time= co2f_io.time.dt.month.isin([10,11]))
co2f_on = co2f_on.resample(time='YE').mean(dim='time')
on_mean = co2f_on.mean(dim='time')

#%% Seasonal climatology of fCO₂ in Indian Ocean 

fig, axs = plt.subplots(2, 2, figsize=(18,14),dpi =600,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace':0.06,'hspace': -0.14})
season = ['DJF','MAM','JJAS',"ON"]
co2f_season = [djf_mean,mam_mean,jjas_mean,on_mean]
xticks = np.arange(35,125,10)
yticks = np.arange(-20,40,10)

from matplotlib import colors
norm = colors.TwoSlopeNorm(vmin=-60,vcenter=0,vmax=60)
season_lv = np.arange(-60,65,5)
slv = np.arange(-60,65,10)
#plotting

for i, ax in enumerate(axs.flat):
    im = ax.contourf(lon,lat,co2f_season[i],levels =season_lv,
                     cmap= cmaps.cmocean_balance,norm=norm,
                     transform=ccrs.PlateCarree(),extend = 'both') 
    ax.coastlines(zorder =15)
    ax.coastlines(zorder = 12)
    ax.set_extent([30, 119, -29, 29], crs=ccrs.PlateCarree())
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
 
    # contours = ax.contour(lon,lat,co2f_season[i],levels = season_lv,
    #                           colors='black', linewidths=0.6)
    # ax.clabel(contours, inline=True, fontsize=12, fmt="%.0f")        

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
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical',ticks =slv )
cbar.ax.tick_params(axis='both', which='major',
               direction='out',   # ticks pointing outward
               length=6,          # tick length
               width=1,labelsize=24)


#titles
cbar.set_label('flux (gC m$^{-2}$ yr$^{-1}$)',fontsize = 24,labelpad=12,
               fontweight = 'bold')
# fig.suptitle(
#     'Seasonal Air–Sea CO$_2$ Flux in Indian Ocean (1994–2024)',
#     fontsize=24,
#     fontweight='bold',
#     y=0.88
# )
plt.show()
#%% standard deviation of seasonal plots 

djf_std = co2f_djf.std(dim =('time'))
jjas_std = co2f_jjas.std(dim =('time'))
mam_std = co2f_mam.std(dim =('time'))
on_std = co2f_on.std(dim =('time'))

fig, axs = plt.subplots(2, 2, figsize=(18,14),
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace':0.06,'hspace': -0.09})
season = ['DJF','MAM','JJAS',"ON"]
season_std = [djf_std,mam_std,jjas_std,on_std]

std_lv = np.arange(0,16,1)
stlv = np.arange(0,16,1)
#plotting

for i, ax in enumerate(axs.flat):
    im = ax.contourf(lon,lat,season_std[i],levels =std_lv,
                     cmap= 'turbo',
                     transform=ccrs.PlateCarree(),extend = 'both') 
    ax.coastlines(zorder =15)
    ax.coastlines(zorder = 12)
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
 
    # contours = ax.contour(lon,lat,co2f_season[i],levels = season_lv,
    #                           colors='black', linewidths=0.6)
    # ax.clabel(contours, inline=True, fontsize=12, fmt="%.0f")        

# ticks and gridlines

    ax.set_xticks(xticks, crs=ccrs.PlateCarree()) # for adding ticks 
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    
    ax.tick_params(axis='both',
                   direction='out',   # makes arrow-like ticks
                   length=5,
                   width=1.2,
                   labelsize=20)
    
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
               width=1,labelsize=20)


#titles
cbar.set_label('flux (gC m$^{-2}$ yr$^{-1}$)',fontsize = 20,labelpad=12,
               fontweight = 'bold')
# fig.suptitle(
#     'Seasonal Air–Sea CO$_2$ Flux in Indian Ocean (1994–2024)',
#     fontsize=24,
#     fontweight='bold',
#     y=0.88
# )
plt.show()

#%% standaed deviation of seasonal anomalies fCO2 in indian ocean

djf_ano = co2f_djf - djf_mean
mam_ano = co2f_mam - mam_mean
jjas_ano = co2f_jjas - jjas_mean
on_ano = co2f_on - on_mean

djf_ano_std = djf_ano.std(dim=('time'))
mam_ano_std = mam_ano.std(dim=('time'))
jjas_ano_std = jjas_ano.std(dim=('time'))
on_ano_std = on_ano.std(dim=('time'))


fig, axs = plt.subplots(2, 2, figsize=(18,14),
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace':0.06,'hspace': -0.09})
season = ['DJF','MAM','JJAS',"ON"]
season_ano_std = [djf_ano_std,mam_ano_std,jjas_ano_std,on_ano_std]

std_ano = np.arange(0,16,1)
ano_lv = np.arange(0,16,1)
#plotting

for i, ax in enumerate(axs.flat):
    im = ax.contourf(lon,lat,season_ano_std[i],levels =std_ano,
                     cmap= 'turbo',
                     transform=ccrs.PlateCarree(),extend = 'both') 
    ax.coastlines(zorder =15)
    ax.coastlines(zorder = 12)
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
 
    # contours = ax.contour(lon,lat,co2f_season[i],levels = season_lv,
    #                           colors='black', linewidths=0.6)
    # ax.clabel(contours, inline=True, fontsize=12, fmt="%.0f")        

# ticks and gridlines

    ax.set_xticks(xticks, crs=ccrs.PlateCarree()) # for adding ticks 
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    
    ax.tick_params(axis='both',
                   direction='out',   # makes arrow-like ticks
                   length=5,
                   width=1.2,
                   labelsize=20)
    
      # Hide labels depending on subplot index
    if i in [0, 1]:          # top row
        ax.tick_params(labelbottom=False)
    
    if i in [1, 3]:          # right column
        ax.tick_params(labelleft=False)

#colorbar

cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.70]) 
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical',ticks =ano_lv )
cbar.ax.tick_params(axis='both', which='major',
               direction='out',   # ticks pointing outward
               length=6,          # tick length
               width=1,labelsize=20)


#titles
cbar.set_label('flux (gC m$^{-2}$ yr$^{-1}$)',fontsize = 20,labelpad=12,
               fontweight = 'bold')
# fig.suptitle(
#     'Seasonal Air–Sea CO$_2$ Flux in Indian Ocean (1994–2024)',
#     fontsize=24,
#     fontweight='bold',
#     y=0.88
# )
plt.show()

#%% Standard deviation of Seasonal pCO₂ in Indian Ocean (1980 -2019)

mon_mean = co2f_io.groupby('time.month').mean()

mon_ano = co2f_io.groupby('time.month') - mon_mean

std_mon = mon_ano.groupby('time.month').std('time')

#%% monthly fco2
fig, axs = plt.subplots(3, 4, figsize=(18, 14),
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace': 0.08, 'hspace': -0.55})
x1 =np.arange(35,120,15)
y1 = np.arange(-30,40,15)
months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

month_lv = np.arange(-60, 65, 5)
mlv = np.arange(-60, 65, 10)

for i, ax in enumerate(axs.flat):

    im = ax.contourf(lon, lat, mon_mean[i],
                     levels=month_lv,
                     cmap=cmaps.cmocean_balance_r,
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

cbar.set_label('flux (gC m$^{-2}$ yr$^{-1}$)',
               fontsize=18,
               labelpad=10,
               fontweight='bold')

# overall spacing
fig.subplots_adjust(left=0.06, right=0.90, top=0.95, bottom=0.07)
fig.suptitle(
    'Monthly Air–Sea CO$_2$ Flux in Indian Ocean (1994–2024)',
    fontsize=24,
    fontweight='bold',
    y=0.87
)
plt.show()
#%% std of fco2
fig, axs = plt.subplots(3, 4, figsize=(18, 14),
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        gridspec_kw={'wspace': 0.08, 'hspace': -0.55})
x1 =np.arange(35,120,15)
y1 = np.arange(-30,40,15)
months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

month_lv = np.arange(0,21,1)
mlv = np.arange(0,21,1)

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

cbar.set_label('flux (gC m$^{-2}$ yr$^{-1}$)',
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