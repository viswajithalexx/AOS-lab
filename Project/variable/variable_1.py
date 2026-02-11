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

#%%
file2 = '/home/bobco-08/24cl05012/CO2/data/data_1/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_multi-vars_35.00E-120.00E_30.00S-30.00N_0.51m_1993-01-01-2024-12-01.nc'

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

        
        mam = da.sel(time = da.time.dt.month.isin([3,4,5])).mean(dim=('time'))
        
        jjas = da.sel(time = da.time.dt.month.isin([6,7,8,9])).mean(dim=('time'))
        
        on = da.sel(time = da.time.dt.month.isin([10,11])).mean(dim=('time'))
        
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
        'levels': np.arange(300,460,10),
        'cmap': cmaps.WhiteBlueGreenYellowRed,
        'title': 'surface ocean pCO$_2$ (1994 - 2024)',
        'unit': 'uatm'

    },
    {
        'fields': {
            'DJF': sss_djf,
            'MAM': sss_mam,
            'JJAS': sss_jjas,
            'ON':  sss_on
        },
        'levels': np.arange(30,38,0.5),
        'cmap': cmaps.MPL_PuBu,
        'title': 'SSS (1994 - 2024)',
        'unit': 'psu'
    },
    {
        'fields': {
            'DJF': sst_djf,
            'MAM': sst_mam,
            'JJAS': sst_jjas,
            'ON':  sst_on
        },
        'levels': np.arange(18,30,1),
        'cmap': cmaps.MPL_RdBu_r,
        'title': 'SST (1994 - 2024)',
        'unit': '$^\circ$C'
    },
    {
        'fields': {
            'DJF': tco2_djf,
            'MAM': tco2_mam,
            'JJAS': tco2_jjas,
            'ON':  tco2_on
        },
        'levels': np.arange(1800,2101,20),
        'cmap': cmaps.MPL_PiYG_r,
        'title': 'Dissolved Inorganic Carbon (1994 - 2024)',
        'unit': 'micro mol kg-1'
    },
    {
        'fields': {
            'DJF': talk_djf,
            'MAM': talk_mam,
            'JJAS': talk_jjas,
            'ON':  talk_on
        },
        'levels': np.arange(2100,2401,20),
        'cmap': cmaps.MPL_PuOr,
        'title': 'Alkalinity (1994 - 2024)',
        'unit': 'micro mol kg-1'
    }
]

import matplotlib.patches as patches

regions = {
    'NWIO': {
        'lat': (5, 25),
        'lon': (45, 65),
        'color': 'k'
    },
    'KC': {
        'lat': (6, 10.5),
        'lon': (74.4, 80.5),
        'color': 'k'
    },
    'NBoB': {
        'lat': (18, 22),
        'lon': (84.5, 94.5),
        'color': 'k'
    }
}

# ---- Main plotting loop ----
for cfg in plot_configs:

    fig, axs = plt.subplots(
        2, 2, figsize=(22,18),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': 0.02, 'hspace': -0.14}
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
        #---- Add region boxes ----
        # for name, reg in regions.items():

        #     lon_min, lon_max = reg['lon']
        #     lat_min, lat_max = reg['lat']

        #     rect = patches.Rectangle(
        #         (lon_min, lat_min),
        #         lon_max - lon_min,
        #         lat_max - lat_min,
        #         linewidth=2,
        #         linestyle = '--',
        #         edgecolor=reg['color'],
        #         facecolor='none',
        #         transform=ccrs.PlateCarree(),
        #         zorder=20
        #     )

        #     ax.add_patch(rect)

        #     # Optional label
        #     ax.text(
        #         lon_min + 0.5,
        #         lat_max + 0.5,
        #         name,
        #         color=reg['color'],
        #         fontsize=12,
        #         fontweight='bold',
        #         transform=ccrs.PlateCarree(),
        #         zorder=21
        #     )

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
        gl.xlabel_style = {'fontsize': 24}
        gl.ylabel_style = {'fontsize': 24}

        if i in [0, 1]:
            gl.bottom_labels = False
        if i in [1, 3]:
            gl.left_labels = False

    # ---- Colorbar ----
    cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks = cfg['levels'])
    cbar.set_label(cfg['unit'], fontsize= 24)
    cbar.ax.tick_params(labelsize=24)

    plt.suptitle(
        f'Seasonal climatology of {cfg["title"]}',
        fontsize=21, y=0.87
    )
    
    plt.savefig(
        f'/home/bobco-08/24cl05012/CO2/plot/plots_1/variables/Seasonal_climatology_{cfg["title"].replace(" ", "_")}.png',
        dpi=600,
        bbox_inches='tight'
    )
    
    plt.close()
#%%
def seasonal_standardeviation(da,time_dim = 'time'):
    

        djf =  da.resample(time='QS-DEC').mean()
        
        djf = djf.isel(time=slice(1, None)) 
                
        djf = djf.sel(time= djf.time.dt.month == 12)
                
        djf = djf.std(dim=('time'))

        
        mam = da.sel(time = da.time.dt.month.isin([3,4,5])).std(dim=('time'))
        
        jjas = da.sel(time = da.time.dt.month.isin([6,7,8,9])).std(dim=('time'))
        
        on = da.sel(time = da.time.dt.month.isin([10,11])).std(dim=('time'))
        
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
        'title': 'surface ocean pCO$_2$ (1994 - 2024)',
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
        'cmap': cmaps.MPL_PuBu,
        'title': 'SSS (1994 - 2024)',
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
        'cmap': cmaps.MPL_RdBu_r,
        'title': 'SST (1994 - 2024)',
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
        'cmap': cmaps.MPL_PiYG_r,
        'title': 'Dissolved Inorganic Carbon (1994 - 2024)',
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
        'cmap': cmaps.MPL_PuOr,
        'title': 'Alkalinity (1994 - 2024)',
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
            levels= 15,
            colors='k',
            linewidths=1
        )
        
        ax.clabel(cs, inline=True, fontsize=8, fmt='%d')

        ax.coastlines(zorder=15)
        ax.add_feature(cf.LAND, facecolor='lightgrey', zorder=10)

        ax.text(
            0.75, 0.98, f"STD ({season})",
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
        f'Seasonal standard deviation  of {cfg["title"]}',
        fontsize=21, y=0.94
    )

    plt.savefig(
        f'/home/bobco-08/24cl05012/CO2/plot/plots_1/variables/Seasonal_std_{cfg["title"].replace(" ", "_")}.png',
        dpi=600,
        bbox_inches='tight'
    )
    
    plt.close()
#%%

ds3 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_thetao_oras_35.00E-112.00E_30.00S-30.00N_0.51-108.03m_1993-01-01-2024-12-01.nc')
ds4 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/variables/cmems_mod_glo_bgc-car_anfc_0.25deg_P1M-m_dissic-talk_35.00E-120.00E_30.00S-30.00N_0.49-222.48m_2021-10-01-2025-12-01.nc')
# ds5  = xr.open_dataset('')

#%%
t_nwio = ds3['thetao_oras'].sel(time= slice('1994-01-01','2024-12-01'),latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim=('latitude','longitude'))

dic_nwio = ds4['dissic'].sel(depth = slice(0.51,100),latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim=('latitude','longitude'))
alk_nwio = ds4['talk'].sel(depth = slice(0.51,100),latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim=('latitude','longitude'))

#%%

t_nwio_monthly = (t_nwio.groupby('time.month').mean('time'))
dic_nwio_monthly = (dic_nwio.groupby('time.month').mean('time'))
alk_nwio_monthly = (alk_nwio.groupby('time.month').mean('time'))


#%%

# nwio_dic = seasonal_climatology(dic_nwio, time_dim='time')

# nwio_dic_djf  = nwio_dic['DJF']
# nwio_dic_mam  = nwio_dic['MAM']
# nwio_dic_jjas = nwio_dic['JJAS']
# nwio_dic_on   = nwio_dic['ON']


# #%%
# nwio_alk = seasonal_climatology(alk_nwio, time_dim='time')

# nwio_alk_djf  = nwio_alk['DJF']
# nwio_alk_mam  = nwio_alk['MAM']
# nwio_alk_jjas = nwio_alk['JJAS']
# nwio_alk_on   = nwio_alk['ON']


#%%
plot_configs = [
    {
        'data': t_nwio_monthly,
        'levels': np.linspace(15,30,20),
        'cmap': 'RdBu_r',
        'title': 'Temperature',
        'unit': '°C'
    },
    # {
    #     'data': dic_nwio_monthly,
    #     'levels': 10,
    #     'cmap': cmaps.NCV_bright,
    #     'title': 'DIC',
    #     'unit': 'µmol kg$^{-1}$'
    # },
    # {
    #     'data': alk_nwio_monthly,
    #     'levels':
    #     'cmap': 'RdBu_r',
    #     'title': 'Alkalinity',
    #     'unit': 'µmol kg$^{-1}$'
    # }
]

month_labels = ['Jan','Feb','Mar','Apr','May','Jun',
                'Jul','Aug','Sep','Oct','Nov','Dec']

for cfg in plot_configs:

    data = cfg['data']

    fig, ax = plt.subplots(figsize=(14, 6))

    cf = ax.contourf(
        data.month,
        data.depth,
        data.transpose('depth', 'month'),
        levels=cfg['levels'],
        cmap=cfg['cmap'],
        extend='both'
    )

    cs = ax.contour(
        data.month,
        data.depth,
        data.transpose('depth', 'month'),
        levels=cfg['levels'],
        colors='k',
        linewidths=0.6
    )

    ax.clabel(cs, fmt='%.1f', fontsize=9)

    ax.invert_yaxis()

    ax.set_xticks(np.arange(1, 13))
    ax.set_xticklabels(month_labels, fontsize=12)

    ax.set_xlabel('Month', fontsize=14)
    ax.set_ylabel('Depth (m)', fontsize=14)

    ax.set_title(
        f'NWIO Monthly Climatology: {cfg["title"]} (Depth–Month)',
        fontsize=18
    )

    cbar = fig.colorbar(cf, ax=ax, pad=0.02)
    cbar.set_label(cfg['unit'], fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    plt.tight_layout()
    plt.show()
#%%

nwio_sst = sst.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))

nwio_sss = sss.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))

nwio_pco2 = pco2.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))

nwio_tco2 = tco2.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))

nwio_talk = talk.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))


#%%

import matplotlib.pyplot as plt
import numpy as np

# ------------------------------
# X-axis: Dec → Nov
# ------------------------------
x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May',
                'Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# ------------------------------
# Monthly climatology
# ------------------------------
t_nwio    = nwio_sst.groupby('time.month').mean('time')
sal_nwio  = nwio_sss.groupby('time.month').mean('time')
pco2_nwio = nwio_pco2.groupby('time.month').mean('time')
dic_nwio  = nwio_tco2.groupby('time.month').mean('time')
alk_nwio  = nwio_talk.groupby('time.month').mean('time')

fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# pCO2 (left axis)
ln1 = ax.plot(x, pco2_nwio, lw=2, color='firebrick', label='pCO$_2$')
ax.set_ylabel('pCO$_2$ (µatm)', color='firebrick')
ax.tick_params(axis='y', labelcolor='firebrick')

# DIC (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, dic_nwio, lw=2, color='purple', label='DIC',linestyle='--')
ax2.set_ylabel('DIC (µmol kg$^{-1}$)', color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# ALK (third axis)
ax3 = ax.twinx()
ax3.spines['right'].set_position(('outward', 60))
ln3 = ax3.plot(x, alk_nwio, lw=2, color='black', label='ALK',linestyle='--')
ax3.set_ylabel('Alkalinity (µmol kg$^{-1}$)', color='black')
ax3.tick_params(axis='y', labelcolor='black')

# Legend
lines = ln1 + ln2 + ln3
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='best')

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
ax.set_title('Northwestern Indian Ocean  Monthly Climatology – Carbonate System')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%%

fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# SST (left axis)
ln1 = ax.plot(x, t_nwio, lw=2, color='darkgreen', label='SST')
ax.set_ylabel('SST (°C)', color='darkgreen')
ax.tick_params(axis='y', labelcolor='darkgreen')

# SSS (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sal_nwio, lw=2, linestyle='--',
               color='navy', label='SSS')
ax2.set_ylabel('SSS (psu)', color='navy')
ax2.tick_params(axis='y', labelcolor='navy')

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='best')

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
ax.set_title('Northwestern Indian Ocean Monthly Climatology – Physical Variables')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()



#%%
mu_sst = sst.sel(latitude = slice(6,10.5),longitude = slice(74.4,80.5)).mean(dim= ('latitude','longitude'))

mu_sss = sss.sel(latitude = slice(6,10.5),longitude = slice(74.4,80.5)).mean(dim= ('latitude','longitude'))

mu_pco2 = pco2.sel(latitude = slice(6,10.5),longitude = slice(74.4,80.5)).mean(dim= ('latitude','longitude'))

mu_tco2 = tco2.sel(latitude = slice(6,10.5),longitude = slice(74.4,80.5)).mean(dim= ('latitude','longitude'))

mu_talk = talk.sel(latitude = slice(6,10.5),longitude = slice(74.4,80.5)).mean(dim= ('latitude','longitude'))

#%%

x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May',
                'Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# Month order for Dec–Nov
months_dec = [1,2,3,4,5,6,7,8,9,10,11,12]

sst_kc  = mu_sst.groupby('time.month').mean('time')
sss_kc  = mu_sss.groupby('time.month').mean('time')
pco2_kc = mu_pco2.groupby('time.month').mean('time')
dic_kc  = mu_tco2.groupby('time.month').mean('time')
alk_kc  = mu_talk.groupby('time.month').mean('time')

# ------------------------------
# Plot
# ------------------------------
fig, ax = plt.subplots(figsize=(11,5),dpi = 600)

# pCO2 (left)
ln1 = ax.plot(x, pco2_kc, lw=2, color='firebrick', label='pCO$_2$')
ax.set_ylabel('pCO$_2$ (µatm)', color='firebrick')
ax.tick_params(axis='y', labelcolor='firebrick')

# DIC (right)
ax2 = ax.twinx()
ln2 = ax2.plot(x, dic_kc, lw=2, color='purple', label='DIC',linestyle='--')
ax2.set_ylabel('DIC (µmol kg$^{-1}$)', color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# ALK (third axis)
ax3 = ax.twinx()
ax3.spines['right'].set_position(('outward', 60))
ln3 = ax3.plot(x, alk_kc, lw=2, color='black', label='ALK',linestyle='--')
ax3.set_ylabel('Alkalinity (µmol kg$^{-1}$)', color='black')
ax3.tick_params(axis='y', labelcolor='black')

# Legend
lines = ln1 + ln2 + ln3
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='best')

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
ax.set_title('Southern Region Monthly Climatology – Carbonate System')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%%

fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# SST (left)
ln1 = ax.plot(x, sst_kc, lw=2, color='navy', label='SST')
ax.set_ylabel('SST (°C)', color='navy')
ax.tick_params(axis='y', labelcolor='navy')

# SSS (right)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_kc, lw=2, linestyle='--',
               color='darkgreen', label='SSS')
ax2.set_ylabel('SSS (psu)', color='darkgreen')
ax2.tick_params(axis='y', labelcolor='darkgreen')

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='best')

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
ax.set_title('Southern Region Monthly Climatology – Physical Variables')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()


#%%

nbob_sst = sst.sel(latitude = slice(18,22),longitude = slice(84.5,94.5)).mean(dim= ('latitude','longitude'))

nbob_sss = sss.sel(latitude = slice(18,22),longitude = slice(84.5,94.5)).mean(dim= ('latitude','longitude'))

nbob_pco2 = pco2.sel(latitude = slice(18,22),longitude = slice(84.5,94.5)).mean(dim= ('latitude','longitude'))

nbob_tco2 = tco2.sel(latitude = slice(18,22),longitude = slice(84.5,94.5)).mean(dim= ('latitude','longitude'))

nbob_talk = talk.sel(latitude = slice(18,22),longitude = slice(84.5,94.5)).mean(dim= ('latitude','longitude'))

#%%

x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May',
                'Jun','Jul','Aug','Sep','Oct','Nov','Dec']

sst_nbob  = nbob_sst.groupby('time.month').mean('time')
sss_nbob  = nbob_sss.groupby('time.month').mean('time')
pco2_nbob = nbob_pco2.groupby('time.month').mean('time')
dic_nbob  = nbob_tco2.groupby('time.month').mean('time')
alk_nbob  = nbob_talk.groupby('time.month').mean('time')


fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# pCO2 (left axis)
ln1 = ax.plot(x, pco2_nbob, lw=2, color='firebrick', label='pCO$_2$')
ax.set_ylabel('pCO$_2$ (µatm)', color='firebrick')
ax.tick_params(axis='y', labelcolor='firebrick')

# DIC (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, dic_nbob, lw=2, color='purple', label='DIC',linestyle='--')
ax2.set_ylabel('DIC (µmol kg$^{-1}$)', color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# ALK (third axis)
ax3 = ax.twinx()
ax3.spines['right'].set_position(('outward', 60))
ln3 = ax3.plot(x, alk_nbob, lw=2, color='black', label='ALK',linestyle='--')
ax3.set_ylabel('ALK (µmol kg$^{-1}$)', color='black')
ax3.tick_params(axis='y', labelcolor='black')

# Legend
lines = ln1 + ln2 + ln3
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='best')

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
ax.set_title('NBoB Monthly Climatology – Carbonate Variables')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%%
fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# SST (left axis)
ln1 = ax.plot(x, sst_nbob, lw=2, color='navy', label='SST')
ax.set_ylabel('SST (°C)', color='navy')
ax.tick_params(axis='y', labelcolor='navy')

# SSS (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_nbob, lw=2, linestyle='--',
               color='darkgreen', label='SSS')
ax2.set_ylabel('SSS (psu)', color='darkgreen')
ax2.tick_params(axis='y', labelcolor='darkgreen')

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='best')

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
ax.set_title('NBoB Monthly Climatology – Physical Variables')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()

#%%

# Convert to 1D arrays (important)
x = nbob_tco2.values
y = nbob_pco2.values

# Remove NaNs (very important for correlation)
mask = np.isfinite(x) & np.isfinite(y)
x = x[mask]
y = y[mask]

# Pearson correlation
r = np.corrcoef(x, y)[0, 1]

# ------------------------------
# Scatter plot
# ------------------------------
plt.figure(figsize=(6, 5))
plt.scatter(x, y, s=40, color='firebrick', alpha=0.7)
plt.xlabel('DIC (µmol kg$^{-1}$)')
plt.ylabel('pCO$_2$ (µatm)')
plt.title(f'NWIO: pCO$_2$ vs DIC (r = {r:.2f})')

plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.show()
