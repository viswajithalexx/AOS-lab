#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 15:35:44 2026

@author: bobco-08
"""

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import cmaps
file1 = '/home/bobco-08/24cl05012/CO2/data/data_1/cmems_obs-mob_glo_bgc-car_my_irr-i_multi-vars_35.12E-119.88E_29.88S-29.88N_1994-01-01-2024-12-01.nc'

ds1 = xr.open_dataset(file1)

#coordinates

lat = ds1['latitude']
lon = ds1['longitude']
t = ds1['time']
#variables
pco2 = ds1['spco2']
tco2 = ds1['tco2']
talk = ds1['talk']

#%%  seasonal time average - pco2

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

#%% tco2

# ---- Seasons & data ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']
seasonal_fields = {
    'DJF': tco2_djf,
    'MAM': tco2_mam,
    'JJAS': tco2_jjas,
    'ON':  tco2_on
}
tco2_levels=  np.arange(1800,2101,25)

fig, axs = plt.subplots(
    2, 2, figsize=(18, 12), dpi=600,
    subplot_kw={'projection': ccrs.PlateCarree()},
    gridspec_kw={'wspace': 0.01, 'hspace': 0.1}
)

for i, ax in enumerate(axs.flat):

    season = seasons[i]
    data = seasonal_fields[season]

    im = ax.contourf(
        data.longitude, data.latitude, data,
        levels= tco2_levels,
        cmap='RdBu_r',
        transform=ccrs.PlateCarree(),
        extend='both'
    )

    ax.coastlines(zorder=15)
    ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

    ax.text(
        0.75, 0.98, f"Mean ({season})",
        transform=ax.transAxes,
        fontsize=18, fontweight='bold',
        va='top', ha='right'
    )
    ax.set_title(
        seasons[i],
        fontsize=18,
        fontweight='bold',
        loc='center'
    )

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
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax,ticks = tco2_levels)
cbar.set_label('micro mol kg-1', fontsize=20)
cbar.ax.tick_params(labelsize=15)

plt.suptitle('Seasonal climatology of TCO$_2$', fontsize=21, y=0.96)

# plt.savefig('tco2_seasonal_climatology.png', dpi=600, bbox_inches='tight')
plt.show()
#%% alk

# ---- Seasons & data ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']
seasonal_fields = {
    'DJF': talk_djf,
    'MAM': talk_mam,
    'JJAS':talk_jjas,
    'ON':  talk_on
}
alk_levels=  np.arange(2100,2401,25)

fig, axs = plt.subplots(
    2, 2, figsize=(18, 12), dpi=600,
    subplot_kw={'projection': ccrs.PlateCarree()},
    gridspec_kw={'wspace': 0.01, 'hspace': 0.1}
)

for i, ax in enumerate(axs.flat):

    season = seasons[i]
    data = seasonal_fields[season]

    im = ax.contourf(
        data.longitude, data.latitude, data,
        levels= alk_levels,
        cmap='RdBu_r',
        transform=ccrs.PlateCarree(),
        extend='both'
    )

    ax.coastlines(zorder=15)
    ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

    ax.text(
        0.75, 0.98, f"Mean ({season})",
        transform=ax.transAxes,
        fontsize=18, fontweight='bold',
        va='top', ha='right'
    )
    ax.set_title(
        seasons[i],
        fontsize=18,
        fontweight='bold',
        loc='center'
    )

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
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax,ticks =alk_levels)
cbar.set_label('micro mol kg-1', fontsize=20)
cbar.ax.tick_params(labelsize=15)

plt.suptitle('Seasonal climatology of TALK', fontsize=21, y=0.96)

# plt.savefig('tco2_seasonal_climatology.png', dpi=600, bbox_inches='tight')
plt.show()

#%%

file2 = '/home/bobco-08/24cl05012/CO2/data/data_1/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_multi-vars_35.00E-120.00E_30.00S-30.00N_0.51m_1993-01-01-2024-12-01.nc'

ds2 = xr.open_dataset(file2)



#variables
sst = ds2['thetao_oras'].sel(time= slice('1994-01-01','2024-12-01')).mean(dim=('depth'))
sss = ds2['so_oras'].sel(time= slice('1994-01-01','2024-12-01')).mean(dim=('depth'))

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
# ---- Seasons & data ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']
seasonal_fields = {
    'DJF': sst_djf,
    'MAM': sst_mam,
    'JJAS':sst_jjas,
    'ON':  sst_on
}

sst_levels=  np.arange(18,31,1)

fig, axs = plt.subplots(
    2, 2, figsize=(18, 12),
    subplot_kw={'projection': ccrs.PlateCarree()},
    gridspec_kw={'wspace': 0.01, 'hspace': 0.1}
)

for i, ax in enumerate(axs.flat):

    season = seasons[i]
    data = seasonal_fields[season]

    im = ax.contourf(
        data.longitude, data.latitude, data,
        levels= sst_levels,
        cmap= 'RdBu_r',
        transform=ccrs.PlateCarree(),
        extend='both'
    )

    ax.coastlines(zorder=15)
    ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

    ax.text(
        0.75, 0.98, f"Mean ({season})",
        transform=ax.transAxes,
        fontsize=18, fontweight='bold',
        va='top', ha='right'
    )
    ax.set_title(
        seasons[i],
        fontsize=18,
        fontweight='bold',
        loc='center'
    )

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
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax,ticks = sst_levels)
cbar.set_label('$^\circ$C', fontsize=20)
cbar.ax.tick_params(labelsize=15)

plt.suptitle('Seasonal climatology of SST', fontsize=21, y=0.96)

# plt.savefig('tco2_seasonal_climatology.png', dpi=600, bbox_inches='tight')
plt.show()
#%%

# ---- Seasons & data ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']
seasonal_fields = {
    'DJF': sss_djf,
    'MAM': sss_mam,
    'JJAS':sss_jjas,
    'ON':  sss_on
}

sss_levels=  np.arange(32,36.2,0.3)

fig, axs = plt.subplots(
    2, 2, figsize=(18, 12),
    subplot_kw={'projection': ccrs.PlateCarree()},
    gridspec_kw={'wspace': 0.01, 'hspace': 0.1}
)

for i, ax in enumerate(axs.flat):

    season = seasons[i]
    data = seasonal_fields[season]

    im = ax.contourf(
        data.longitude, data.latitude, data,
        levels= sss_levels,
        cmap= cmaps.MPL_Blues,
        transform=ccrs.PlateCarree(),
        extend='both'
    )

    ax.coastlines(zorder=15)
    ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

    ax.text(
        0.75, 0.98, f"Mean ({season})",
        transform=ax.transAxes,
        fontsize=18, fontweight='bold',
        va='top', ha='right'
    )
    ax.set_title(
        seasons[i],
        fontsize=18,
        fontweight='bold',
        loc='center'
    )

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
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax,ticks = sss_levels )
cbar.set_label('PSU', fontsize=20)
cbar.ax.tick_params(labelsize=15)

plt.suptitle('Seasonal climatology of SSS', fontsize=21, y=0.96)

# plt.savefig('tco2_seasonal_climatology.png', dpi=600, bbox_inches='tight')
plt.show()
#%%

nwio_t =  ds2['thetao_oras'].sel(time= slice('1994-01-01','2024-12-01'), latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim=('depth'))

nwio_sst = seasonal_climatology(nwio_t, time_dim='time')

nwio_djf  = nwio_sst['DJF']
nwio_mam  = nwio_sst['MAM']
nwio_jjas = nwio_sst['JJAS']
nwio_on   = nwio_sst['ON']

#%%
nwio_s =  ds2['so_oras'].sel(time= slice('1994-01-01','2024-12-01'), latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim=('depth'))

nwio_sss = seasonal_climatology(nwio_s, time_dim='time')

nwio_djf_sss  = nwio_sss['DJF']
nwio_mam_sss   = nwio_sss['MAM']
nwio_jjas_sss  = nwio_sss['JJAS']
nwio_on_sss    = nwio_sss['ON']
#%%

nwio_talk =  ds1['talk'].sel(time= slice('1994-01-01','2024-12-01'), latitude = slice(5,22.5),longitude = slice(45,65))

nwio_alk = seasonal_climatology(nwio_talk, time_dim='time')

nwio_djf_alk  = nwio_alk['DJF']
nwio_mam_alk   = nwio_alk['MAM']
nwio_jjas_alk  = nwio_alk['JJAS']
nwio_on_alk    = nwio_alk['ON']
#%%

nwio_dic=  ds1['tco2'].sel(time= slice('1994-01-01','2024-12-01'), latitude = slice(5,22.5),longitude = slice(45,65))

nwio_tco2 = seasonal_climatology(nwio_dic, time_dim='time')

nwio_djf_tco2  = nwio_tco2['DJF']
nwio_mam_tco2   = nwio_tco2['MAM']
nwio_jjas_tco2  = nwio_tco2['JJAS']
nwio_on_tco2    = nwio_tco2['ON']

#%%
seasons = ['DJF', 'MAM', 'JJAS', 'ON']
#%%
plot_configs = [
    {
        'seasonal_fields': {
            'DJF': nwio_djf_tco2,
            'MAM': nwio_mam_tco2,
            'JJAS': nwio_jjas_tco2,
            'ON':  nwio_on_tco2
        },
        'levels': np.linspace(1980, 2090, 15),
        'cl': cmaps.rainbow,
        'title': 'DIC',
        'unit': 'micro mol kg-1'
    },
    {
        'seasonal_fields': {
            'DJF': nwio_djf_alk,
            'MAM': nwio_mam_alk,
            'JJAS': nwio_jjas_alk,
            'ON':  nwio_on_alk
        },
        'levels': np.linspace(2320, 2390, 15),
        'cl': cmaps.rainbow,
        'title': 'ALK',
        'unit': 'micro mol kg-1'
    },
    {
        'seasonal_fields': {
            'DJF': nwio_djf,
            'MAM': nwio_mam,
            'JJAS': nwio_jjas,
            'ON':  nwio_on
        },
        'levels': np.linspace(22, 30, 15),
        'cl': cmaps.cmp_b2r,
        'title': 'SST',
        'unit': '$^\circ$C'
    },
    {
        'seasonal_fields': {
            'DJF': nwio_djf_sss,
            'MAM': nwio_mam_sss,
            'JJAS': nwio_jjas_sss,
            'ON':  nwio_on_sss
        },
        'levels': np.round(np.linspace(34, 37, 15), 2),
        'cl': cmaps.rainbow,
        'title': 'SSS',
        'unit': 'PSU'
    }
]

for cfg in plot_configs:

    seasonal_fields = cfg['seasonal_fields']
    levels = cfg['levels']
    cl = cfg['cl']
    title = cfg['title']
    unit = cfg['unit']

    fig, axs = plt.subplots(
        2, 2, figsize=(20, 18),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': 0.05, 'hspace': 0.1}
    )

    for i, ax in enumerate(axs.flat):

        season = seasons[i]
        data = seasonal_fields[season]

        im = ax.contourf(
            data.longitude, data.latitude, data,
            levels=levels,
            cmap=cl,
            transform=ccrs.PlateCarree(),
            extend='both'
        )

        cs = ax.contour(
            data.longitude, data.latitude, data,
            levels=levels,
            colors='k',
            linewidths=0.6,
            transform=ccrs.PlateCarree()
        )

        ax.clabel(cs, fmt="%.1f", fontsize=8)

        ax.coastlines(zorder=15)
        ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

        ax.text(
            0.75, 0.98, f"Mean ({season})",
            transform=ax.transAxes,
            fontsize=18, fontweight='bold',
            va='top', ha='right'
        )

        ax.set_title(season, fontsize=18, fontweight='bold')

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
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=levels)
    cbar.set_label(unit, fontsize=20)
    cbar.ax.tick_params(labelsize=15)

    plt.suptitle(f'Seasonal climatology of NWIO {title}', fontsize=21, y=0.94)

    plt.show()

#%%
nas_t = ds2['thetao_oras'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(22.5, 28),
    longitude=slice(56, 70)
).mean(dim=('depth'))

nas_sst = seasonal_climatology(nas_t, time_dim='time')

nas_djf  = nas_sst['DJF']
nas_mam  = nas_sst['MAM']
nas_jjas = nas_sst['JJAS']
nas_on   = nas_sst['ON']

#%%
nas_s = ds2['so_oras'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(22.5, 28),
    longitude=slice(56, 70)
).mean(dim=('depth'))

nas_sss = seasonal_climatology(nas_s, time_dim='time')

nas_djf_sss   = nas_sss['DJF']
nas_mam_sss   = nas_sss['MAM']
nas_jjas_sss  = nas_sss['JJAS']
nas_on_sss    = nas_sss['ON']

#%%
nas_talk = ds1['talk'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(22.5, 28),
    longitude=slice(56, 70)
)

nas_alk = seasonal_climatology(nas_talk, time_dim='time')

nas_djf_alk   = nas_alk['DJF']
nas_mam_alk   = nas_alk['MAM']
nas_jjas_alk  = nas_alk['JJAS']
nas_on_alk    = nas_alk['ON']

#%%
nas_dic = ds1['tco2'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(22.5, 28),
    longitude=slice(56, 70)
)

nas_tco2 = seasonal_climatology(nas_dic, time_dim='time')

nas_djf_tco2   = nas_tco2['DJF']
nas_mam_tco2   = nas_tco2['MAM']
nas_jjas_tco2  = nas_tco2['JJAS']
nas_on_tco2    = nas_tco2['ON']

#%%
seasons = ['DJF', 'MAM', 'JJAS', 'ON']
#%%
plot_configs = [
    {
        'seasonal_fields': {
            'DJF': nas_djf_tco2,
            'MAM': nas_mam_tco2,
            'JJAS': nas_jjas_tco2,
            'ON':  nas_on_tco2
        },
        'levels': np.linspace(1980, 2090, 15),
        'cl': cmaps.rainbow,
        'title': 'DIC',
        'unit': 'micro mol kg-1'
    },
    {
        'seasonal_fields': {
            'DJF': nas_djf_alk,
            'MAM': nas_mam_alk,
            'JJAS': nas_jjas_alk,
            'ON':  nas_on_alk
        },
        'levels': np.linspace(2320, 2390, 15),
        'cl': cmaps.rainbow,
        'title': 'ALK',
        'unit': 'micro mol kg-1'
    },
    {
        'seasonal_fields': {
            'DJF': nas_djf,
            'MAM': nas_mam,
            'JJAS': nas_jjas,
            'ON':  nas_on
        },
        'levels': np.linspace(22, 30, 15),
        'cl': cmaps.cmp_b2r,
        'title': 'SST',
        'unit': '$^\circ$C'
    },
    {
        'seasonal_fields': {
            'DJF': nas_djf_sss,
            'MAM': nas_mam_sss,
            'JJAS': nas_jjas_sss,
            'ON':  nas_on_sss
        },
        'levels': np.round(np.linspace(34, 37, 15), 2),
        'cl': cmaps.rainbow,
        'title': 'SSS',
        'unit': 'PSU'
    }
]


for cfg in plot_configs:

    seasonal_fields = cfg['seasonal_fields']
    levels = cfg['levels']
    cl = cfg['cl']
    title = cfg['title']
    unit = cfg['unit']

    fig, axs = plt.subplots(
        2, 2, figsize=(15,12),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': 0.05, 'hspace': -0.6}
    )

    for i, ax in enumerate(axs.flat):

        season = seasons[i]
        data = seasonal_fields[season]

        im = ax.contourf(
            data.longitude, data.latitude, data,
            levels=levels,
            cmap=cl,
            transform=ccrs.PlateCarree(),
            extend='both'
        )

        cs = ax.contour(
            data.longitude, data.latitude, data,
            levels=levels,
            colors='k',
            linewidths=0.6,
            transform=ccrs.PlateCarree()
        )

        ax.clabel(cs, fmt="%.1f", fontsize=8)

        ax.coastlines(zorder=15)
        ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

        ax.text(
            0.75, 0.98, f"Mean ({season})",
            transform=ax.transAxes,
            fontsize=18, fontweight='bold',
            va='top', ha='right'
        )

        ax.set_title(season, fontsize=18, fontweight='bold')

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
    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.50])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=levels)
    cbar.set_label(unit, fontsize=20)
    cbar.ax.tick_params(labelsize=15)

    plt.suptitle(f'Seasonal climatology of NAS {title}', fontsize=21, y=0.79)

    plt.show()
#%%

eio_t = ds2['thetao_oras'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(-6.5, 5),
    longitude=slice(49, 92)
).mean(dim=('depth'))

eio_sst = seasonal_climatology(eio_t, time_dim='time')

eio_djf  = eio_sst['DJF']
eio_mam  = eio_sst['MAM']
eio_jjas = eio_sst['JJAS']
eio_on   = eio_sst['ON']

#%%
eio_s = ds2['so_oras'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(-6.5, 5),
    longitude=slice(49, 92)
).mean(dim=('depth'))

eio_sss = seasonal_climatology(eio_s, time_dim='time')

eio_djf_sss   = eio_sss['DJF']
eio_mam_sss   = eio_sss['MAM']
eio_jjas_sss  = eio_sss['JJAS']
eio_on_sss    = eio_sss['ON']

#%%
eio_talk = ds1['talk'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(-6.5, 5),
    longitude=slice(49, 92)
)

eio_alk = seasonal_climatology(eio_talk, time_dim='time')

eio_djf_alk   = eio_alk['DJF']
eio_mam_alk   = eio_alk['MAM']
eio_jjas_alk  = eio_alk['JJAS']
eio_on_alk    = eio_alk['ON']

#%%
eio_dic = ds1['tco2'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(-6.5, 5),
    longitude=slice(49, 92)
)

eio_tco2 = seasonal_climatology(eio_dic, time_dim='time')

eio_djf_tco2   = eio_tco2['DJF']
eio_mam_tco2   = eio_tco2['MAM']
eio_jjas_tco2  = eio_tco2['JJAS']
eio_on_tco2    = eio_tco2['ON']

#%%
seasons = ['DJF', 'MAM', 'JJAS', 'ON']
#%%
plot_configs = [
    {
        'seasonal_fields': {
            'DJF': eio_djf_tco2,
            'MAM': eio_mam_tco2,
            'JJAS': eio_jjas_tco2,
            'ON':  eio_on_tco2
        },
        'levels': np.linspace(1980, 2090, 15),
        'cl': cmaps.rainbow,
        'title': 'DIC',
        'unit': 'micro mol kg-1'
    },
    {
        'seasonal_fields': {
            'DJF': eio_djf_alk,
            'MAM': eio_mam_alk,
            'JJAS': eio_jjas_alk,
            'ON':  eio_on_alk
        },
        'levels': np.linspace(2320, 2390, 15),
        'cl': cmaps.rainbow,
        'title': 'ALK',
        'unit': 'micro mol kg-1'
    },
    {
        'seasonal_fields': {
            'DJF': eio_djf,
            'MAM': eio_mam,
            'JJAS': eio_jjas,
            'ON':  eio_on
        },
        'levels': np.linspace(22, 30, 15),
        'cl': cmaps.cmp_b2r,
        'title': 'SST',
        'unit': '$^\circ$C'
    },
    {
        'seasonal_fields': {
            'DJF': eio_djf_sss,
            'MAM': eio_mam_sss,
            'JJAS': eio_jjas_sss,
            'ON':  eio_on_sss
        },
        'levels': np.round(np.linspace(34, 37, 15), 2),
        'cl': cmaps.rainbow,
        'title': 'SSS',
        'unit': 'PSU'
    }
]


for cfg in plot_configs:

    seasonal_fields = cfg['seasonal_fields']
    levels = cfg['levels']
    cl = cfg['cl']
    title = cfg['title']
    unit = cfg['unit']

    fig, axs = plt.subplots(
        2, 2, figsize=(15,12),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': 0.05, 'hspace': -0.7}
    )

    for i, ax in enumerate(axs.flat):

        season = seasons[i]
        data = seasonal_fields[season]

        im = ax.contourf(
            data.longitude, data.latitude, data,
            levels=levels,
            cmap=cl,
            transform=ccrs.PlateCarree(),
            extend='both'
        )

        cs = ax.contour(
            data.longitude, data.latitude, data,
            levels=levels,
            colors='k',
            linewidths=0.6,
            transform=ccrs.PlateCarree()
        )

        ax.clabel(cs, fmt="%.1f", fontsize=8)

        ax.coastlines(zorder=15)
        ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)
    

        ax.set_title(season, fontsize=18, fontweight='bold')

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
    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.50])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=levels)
    cbar.set_label(unit, fontsize=20)
    cbar.ax.tick_params(labelsize=15)

    plt.suptitle(f'Seasonal climatology of EIO {title}', fontsize=21, y=0.72)

    plt.show()
#%%
esio_t = ds2['thetao_oras'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(-6.6, 8),
    longitude=slice(92, 109)
).mean(dim=('depth'))

esio_sst = seasonal_climatology(esio_t, time_dim='time')

esio_djf  = esio_sst['DJF']
esio_mam  = esio_sst['MAM']
esio_jjas = esio_sst['JJAS']
esio_on   = esio_sst['ON']

#%%
esio_s = ds2['so_oras'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(-6.6, 8),
    longitude=slice(92, 109)
).mean(dim=('depth'))

esio_sss = seasonal_climatology(esio_s, time_dim='time')

esio_djf_sss   = esio_sss['DJF']
esio_mam_sss   = esio_sss['MAM']
esio_jjas_sss  = esio_sss['JJAS']
esio_on_sss    = esio_sss['ON']

#%%
esio_talk = ds1['talk'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(-6.6, 8),
    longitude=slice(92, 109)
)

esio_alk = seasonal_climatology(esio_talk, time_dim='time')

esio_djf_alk   = esio_alk['DJF']
esio_mam_alk   = esio_alk['MAM']
esio_jjas_alk  = esio_alk['JJAS']
esio_on_alk    = esio_alk['ON']

#%%
esio_dic = ds1['tco2'].sel(
    time=slice('1994-01-01', '2024-12-01'),
    latitude=slice(-6.6, 8),
    longitude=slice(92, 109)
)

esio_tco2 = seasonal_climatology(esio_dic, time_dim='time')

esio_djf_tco2   = esio_tco2['DJF']
esio_mam_tco2   = esio_tco2['MAM']
esio_jjas_tco2  = esio_tco2['JJAS']
esio_on_tco2    = esio_tco2['ON']

#%%
seasons = ['DJF', 'MAM', 'JJAS', 'ON']

#%%
plot_configs = [
    {
        'seasonal_fields': {
            'DJF': esio_djf_tco2,
            'MAM': esio_mam_tco2,
            'JJAS': esio_jjas_tco2,
            'ON':  esio_on_tco2
        },
        'levels': np.linspace(1834,1930, 15),
        'cl': cmaps.rainbow,
        'title': 'DIC',
        'unit': 'micro mol kg-1'
    },
    {
        'seasonal_fields': {
            'DJF': esio_djf_alk,
            'MAM': esio_mam_alk,
            'JJAS': esio_jjas_alk,
            'ON':  esio_on_alk
        },
        'levels': np.linspace(2125, 2250, 15),
        'cl': cmaps.rainbow,
        'title': 'ALK',
        'unit': 'micro mol kg-1'
    },
    {
        'seasonal_fields': {
            'DJF': esio_djf,
            'MAM': esio_mam,
            'JJAS': esio_jjas,
            'ON':  esio_on
        },
        'levels': np.linspace(25, 30, 15),
        'cl': cmaps.cmp_b2r,
        'title': 'SST',
        'unit': '$^\circ$C'
    },
    {
        'seasonal_fields': {
            'DJF': esio_djf_sss,
            'MAM': esio_mam_sss,
            'JJAS': esio_jjas_sss,
            'ON':  esio_on_sss
        },
        'levels': np.round(np.linspace(30, 34, 15), 2),
        'cl': cmaps.MPL_PuBu,
        'title': 'SSS',
        'unit': 'PSU'
    }
]

for cfg in plot_configs:

    seasonal_fields = cfg['seasonal_fields']
    levels = cfg['levels']
    cl = cfg['cl']
    title = cfg['title']
    unit = cfg['unit']

    fig, axs = plt.subplots(
        2, 2, figsize=(18,15),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': -0.02, 'hspace': 0.2}
    )

    for i, ax in enumerate(axs.flat):

        season = seasons[i]
        data = seasonal_fields[season]

        im = ax.contourf(
            data.longitude, data.latitude, data,
            levels=levels,
            cmap=cl,
            transform=ccrs.PlateCarree(),
            extend='both'
        )

        cs = ax.contour(
            data.longitude, data.latitude, data,
            levels=levels,
            colors='k',
            linewidths=0.6,
            transform=ccrs.PlateCarree()
        )

        ax.clabel(cs, fmt="%.1f", fontsize=8)

        ax.coastlines(zorder=15)
        ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

        ax.set_title(season, fontsize=18, fontweight='bold')

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
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.70])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=levels)
    cbar.set_label(unit, fontsize=20)
    cbar.ax.tick_params(labelsize=15)

    plt.suptitle(f'Seasonal climatology of ESIO {title}', fontsize=21, y=0.95)

    plt.show()
