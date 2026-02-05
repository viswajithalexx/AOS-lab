#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 18:41:37 2026

@author: bobco-08
"""
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import cmaps
file1 = '/home/bobco-08/24cl05012/CO2/data/data_1/decomposition_data/cmems_obs-mob_glo_bgc-car_my_irr-i_multi-vars_35.12E-119.88E_29.88S-29.88N_1994-01-01-2024-12-01.nc'

ds1 = xr.open_dataset(file1)

#coordinates

lat = ds1['latitude']
lon = ds1['longitude']
t = ds1['time']
#variables
pco2 = ds1['spco2']
tco2 = ds1['tco2']
talk = ds1['talk']

#%%
file2 = '/home/bobco-08/24cl05012/CO2/data/data_1/decomposition_data/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_multi-vars_35.00E-120.00E_30.00S-30.00N_0.51m_1993-01-01-2024-12-01.nc'

ds2 = xr.open_dataset(file2)

#variables
sst = ds2['thetao_oras'].sel(time= slice('1994-01-01','2024-12-01')).mean(dim=('depth'))
sss = ds2['so_oras'].sel(time= slice('1994-01-01','2024-12-01')).mean(dim=('depth'))
#%%

def seasonal_climatology(da,time_dim = 'time'):
    

        djf =  da.resample(time='QS-DEC').mean(dim=('time'))
        
        djf = djf.isel(time=slice(1, None)) 
                
        djf = djf.sel(time= djf.time.dt.month == 12)
                
        djf = djf.mean(dim=('time'))

        
        mam = da.sel(time = pco2.time.dt.month.isin([3,4,5])).mean(dim=('time'))
        
        jjas = da.sel(time = pco2.time.dt.month.isin([6,7,8,9])).mean(dim=('time'))
        
        on = da.sel(time = pco2.time.dt.month.isin([10,11])).mean(dim=('time'))
        
        return { 'DJF': djf,
                  'MAM': mam,
                 'JJAS': jjas,
                  'ON': on}    
    
#%%

seasons_pco2 = seasonal_climatology(pco2, time_dim='time')

sp_djf  = seasons_pco2['DJF']
sp_mam  = seasons_pco2['MAM']
sp_jjas = seasons_pco2['JJAS']
sp_on   = seasons_pco2['ON']
   
#%%

seasonal_talk = seasonal_climatology(talk, time_dim='time')

talk_djf  = seasonal_talk['DJF']
talk_mam  = seasonal_talk['MAM']
talk_jjas = seasonal_talk['JJAS']
talk_on   = seasonal_talk['ON']
#%%
seasonal_tco2 = seasonal_climatology(tco2, time_dim='time')

tco2_djf  = seasonal_tco2['DJF']
tco2_mam  = seasonal_tco2['MAM']
tco2_jjas = seasonal_tco2['JJAS']
tco2_on   = seasonal_tco2['ON']

#%%

seasons_sst = seasonal_climatology(sst, time_dim='time')

sst_djf  = seasons_sst['DJF']
sst_mam  = seasons_sst['MAM']
sst_jjas = seasons_sst['JJAS']
sst_on   = seasons_sst['ON']

#%%

seasons_sss = seasonal_climatology(sss, time_dim='time')

sss_djf  = seasons_sss['DJF']
sss_mam  = seasons_sss['MAM']
sss_jjas = seasons_sss['JJAS']
sss_on   = seasons_sss['ON']

#%%
# ---- Seasons ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']

# ---- Plot configurations ----
plot_configs = [
    {
        'fields': {
            'DJF': sp_djf,
            'MAM':  sp_mam,
            'JJAS':  sp_jjas,
            'ON':   sp_on
        },
        'levels': np.linspace(300,480,16,dtype = int),
        'cmap': cmaps.WhiteBlueGreenYellowRed,
        'title': 'surface ocean pCO2',
        'unit': '$^\circ$C'
    },
    
    {
        'fields': {
            'DJF': sst_djf,
            'MAM': sst_mam,
            'JJAS': sst_jjas,
            'ON':  sst_on
        },
        'levels': np.linspace(18,30,15),
        'cmap': cmaps.cmp_b2r,
        'title': 'SSS',
        'unit': 'PSU'
    },
    {
        'fields': {
            'DJF': tco2_djf,
            'MAM': tco2_mam,
            'JJAS': tco2_jjas,
            'ON':  tco2_on
        },
        'levels': np.linspace(1800, 2101, 15),
        'cmap': cmaps.NCV_bright,
        'title': 'TCO$_2$',
        'unit': 'micro mol kg-1'
    },
    {
        'fields': {
            'DJF': talk_djf,
            'MAM': talk_mam,
            'JJAS': talk_jjas,
            'ON':  talk_on
        },
        'levels': np.arange(2100, 2401,20),
        'cmap': 'RdBu_r',
        'title': 'TALK',
        'unit': 'micro mol kg-1'
    }
]

# ---- Main plotting loop ----
for cfg in plot_configs:

    fig, axs = plt.subplots(
        2, 2, figsize=(18, 12),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': -0.10, 'hspace': 0.04}
    )

    for i, ax in enumerate(axs.flat):

        season = seasons[i]
        data = cfg['fields'][season]

        im = ax.contourf(
            data.longitude, data.latitude, data,
            levels=cfg['levels'],
            cmap=cfg['cmap'],
            transform=ccrs.PlateCarree(),
            extend='both'
        )
        
        cs = ax.contour(
            data.longitude, data.latitude, data,
            levels=cfg['levels'],
            colors='k',
            linewidths=1
        )
        
        ax.clabel(cs, inline=True, fontsize=8, fmt='%d')

        ax.coastlines(zorder=15)
        ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

        ax.text(
            0.75, 0.98, f"Mean ({season})",
            transform=ax.transAxes,
            fontsize=18, zorder = 15,fontweight='bold',
            va='top', ha='right'
        )

        # ax.set_title(season, fontsize=18, fontweight='bold')

        gl = ax.gridlines(draw_labels=True, visible=False,
                          linewidth=0.5, color='grey', alpha=0.2)
        gl.right_labels = False
        gl.top_labels = False
        gl.xlabel_style = {'fontsize': 15}
        gl.ylabel_style = {'fontsize': 15}

        if i in [0, 1]:
            gl.bottom_labels = False
        if i in [1, 3]:
            gl.left_labels = False

    # ---- Colorbar ----
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=cfg['levels'])
    cbar.set_label(cfg['unit'], fontsize=20)
    cbar.ax.tick_params(labelsize=15)

    plt.suptitle(
        f'Seasonal climatology of {cfg["title"]}',
        fontsize=21, y=0.94
    )

    plt.show()
#%%
def seasonal_standardeviation(da,time_dim = 'time'):
    

        djf =  sst.resample(time='QS-DEC').mean()
        
        djf = djf.isel(time=slice(1, None)) 
                
        djf = djf.sel(time= djf.time.dt.month == 12)
                
        djf = djf.std(dim=('time'))

        
        mam = da.sel(time = pco2.time.dt.month.isin([3,4,5])).std(dim=('time'))
        
        jjas = da.sel(time = pco2.time.dt.month.isin([6,7,8,9])).std(dim=('time'))
        
        on = da.sel(time = pco2.time.dt.month.isin([10,11])).std(dim=('time'))
        
        return { 'DJF': djf,
                  'MAM': mam,
                 'JJAS': jjas,
                  'ON': on}    
    
#%%

seasons_pco2 = seasonal_standardeviation(pco2, time_dim='time')

sp_djf_std  = seasons_pco2['DJF']
sp_mam_std  = seasons_pco2['MAM']
sp_jjas_std = seasons_pco2['JJAS']
sp_on_std   = seasons_pco2['ON']
   
#%%

seasonal_talk = seasonal_standardeviation(talk, time_dim='time')

talk_djf_std  = seasonal_talk['DJF']
talk_mam_std  = seasonal_talk['MAM']
talk_jjas_std = seasonal_talk['JJAS']
talk_on_std   = seasonal_talk['ON']
#%%
seasonal_tco2 = seasonal_standardeviation(tco2, time_dim='time')

tco2_djf_std  = seasonal_tco2['DJF']
tco2_mam_std  = seasonal_tco2['MAM']
tco2_jjas_std = seasonal_tco2['JJAS']
tco2_on_std   = seasonal_tco2['ON']

#%%

seasons_sst = seasonal_standardeviation(sst, time_dim='time')

sst_djf_std  = seasons_sst['DJF']
sst_mam_std  = seasons_sst['MAM']
sst_jjas_std = seasons_sst['JJAS']
sst_on_std   = seasons_sst['ON']

#%%

seasons_sss = seasonal_standardeviation(sss, time_dim='time')

sss_djf_std  = seasons_sss['DJF']
sss_mam_std  = seasons_sss['MAM']
sss_jjas_std = seasons_sss['JJAS']
sss_on_std   = seasons_sss['ON']

#%%
# ---- Seasons ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']

# ---- Plot configurations ----
plot_configs = [
    {
        'fields': {
            'DJF': sp_djf_std,
            'MAM':  sp_mam_std,
            'JJAS':  sp_jjas_std,
            'ON':   sp_on_std
        },
        'levels': 20,
        'cmap': cmaps.WhiteBlueGreenYellowRed,
        'title': 'surface ocean pCO2',
        'unit': 'uatm'
    },
    {
        'fields': {
            'DJF': sss_djf_std,
            'MAM': sss_mam_std,
            'JJAS': sss_jjas_std,
            'ON':  sss_on_std
        },
        'levels': 20,
        'cmap': cmaps.MPL_Blues,
        'title': 'SSS',
        'unit': 'PSU'
    },
    {
        'fields': {
            'DJF': sst_djf_std,
            'MAM': sst_mam_std,
            'JJAS': sst_jjas_std,
            'ON':  sst_on_std
        },
        'levels': 20,
        'cmap': cmaps.cmp_b2r,
        'title': 'SST',
        'unit': '$^\circ$C'
    },
    {
        'fields': {
            'DJF': tco2_djf_std,
            'MAM': tco2_mam_std,
            'JJAS': tco2_jjas_std,
            'ON':  tco2_on_std
        },
        'levels': 20,
        'cmap': cmaps.NCV_bright,
        'title': 'TCO$_2$',
        'unit': 'micro mol kg-1'
    },
    {
        'fields': {
            'DJF': talk_djf_std,
            'MAM': talk_mam_std,
            'JJAS': talk_jjas_std,
            'ON':  talk_on_std
        },
        'levels': 20,
        'cmap': 'RdBu_r',
        'title': 'TALK',
        'unit': 'micro mol kg-1'
    }
]

# ---- Main plotting loop ----
for cfg in plot_configs:

    fig, axs = plt.subplots(
        2, 2, figsize=(18, 12),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': -0.10, 'hspace': 0.04}
    )

    for i, ax in enumerate(axs.flat):

        season = seasons[i]
        data = cfg['fields'][season]

        im = ax.contourf(
            data.longitude, data.latitude, data,
            levels=cfg['levels'],
            cmap=cfg['cmap'],
            transform=ccrs.PlateCarree(),
            extend='both'
        )
        
        # cs = ax.contour(
        #     data.longitude, data.latitude, data,
        #     levels=cfg['levels'],
        #     colors='k',
        #     linewidths=1
        # )
        
        # ax.clabel(cs, inline=True, fontsize=8, fmt='%d')

        ax.coastlines(zorder=15)
        ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

        ax.text(
            0.75, 0.98, f"Mean ({season})",
            transform=ax.transAxes,
            fontsize=18, zorder = 15,fontweight='bold',
            va='top', ha='right'
        )

        # ax.set_title(season, fontsize=18, fontweight='bold')

        gl = ax.gridlines(draw_labels=True, visible=False,
                          linewidth=0.5, color='grey', alpha=0.2)
        gl.right_labels = False
        gl.top_labels = False
        gl.xlabel_style = {'fontsize': 15}
        gl.ylabel_style = {'fontsize': 15}

        if i in [0, 1]:
            gl.bottom_labels = False
        if i in [1, 3]:
            gl.left_labels = False

    # ---- Colorbar ----
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax,)
    cbar.set_label(cfg['unit'], fontsize=20)
    cbar.ax.tick_params(labelsize=15)

    plt.suptitle(
        f'Seasonal climatology of {cfg["title"]}',
        fontsize=21, y=0.94
    )

    plt.show()
#%%
ds3 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_thetao_oras_35.00E-112.00E_30.00S-30.00N_0.51-108.03m_1993-01-01-2024-12-01.nc')

temp = ds3['thetao_oras'].sel(latitude = slice(5,22.5), longitude = slice(45,65)).mean(dim=('latitude','time'))

lon =    temp.longitude     # longitude
depth =    temp.depth    # depth (m)

# ---- Plot ----
fig, ax = plt.subplots(figsize=(10, 5))

levels = np.linspace(5,27,30)

# Filled contours
cf = ax.contourf(
    lon,
    depth,
    temp,
    cmap = cmaps.WhiteBlueGreenYellowRed,
    levels=levels,
    extend='both'
)

# Contour lines (overlay)
cs = ax.contour(
    lon,
    depth,
    temp,
    levels=levels,
    colors='k',
    linewidths=0.6
)

# Add contour labels
ax.clabel(cs, fmt='%d', fontsize=8)

# Invert y-axis so depth increases downward
ax.invert_yaxis()

# Labels and title
ax.set_xlabel('Longitude (°E)')
ax.set_ylabel('Depth (m)')
ax.set_title('Temperature (°C): Depth vs Longitude')

# Colorbar
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label('Temperature (°C)')

plt.tight_layout()
plt.show()
