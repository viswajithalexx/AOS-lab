#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 18:41:37 2026

@author: bobco-08
"""
#importing library

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import cmaps


file1 = '/home/bobco-08/24cl05012/CO2/data/data_1/cmems/bio_var/cmems_obs-mob_glo_bgc-car_my_irr-i_1775027250539_30E-119.88E_29.88S-29.88N_1994-01-01-2024-12-01.nc'

ds1 = xr.open_dataset(file1)

#coordinates

lat = ds1['latitude']
lon = ds1['longitude']
t = ds1['time']
#variables
pco2 = ds1['spco2']
tco2 = ds1['tco2']
talk = ds1['talk']

file2 = '/home/bobco-08/24cl05012/CO2/data/data_1/cmems/phy_var/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_1775039295106.nc'

ds2 = xr.open_dataset(file2)

#variables
sst = ds2['thetao_oras'].sel(time= slice('1994-01-01','2024-12-01')).mean(dim=('depth'))
sss = ds2['so_oras'].sel(time= slice('1994-01-01','2024-12-01')).mean(dim=('depth'))
#%% seasonal climatology function

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
#%% standard deviation function
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

    
#%% seasonal pco2,sst,sss,alk,dic climatology

seasons_pco2 = seasonal_climatology(pco2, time_dim='time')

sp_djf  = seasons_pco2['DJF']
sp_mam  = seasons_pco2['MAM']
sp_jjas = seasons_pco2['JJAS']
sp_on   = seasons_pco2['ON']
   
seasonal_talk = seasonal_climatology(talk, time_dim='time')

talk_djf  = seasonal_talk['DJF']
talk_mam  = seasonal_talk['MAM']
talk_jjas = seasonal_talk['JJAS']
talk_on   = seasonal_talk['ON']

seasonal_tco2 = seasonal_climatology(tco2, time_dim='time')

tco2_djf  = seasonal_tco2['DJF']
tco2_mam  = seasonal_tco2['MAM']
tco2_jjas = seasonal_tco2['JJAS']
tco2_on   = seasonal_tco2['ON']



seasons_sst = seasonal_climatology(sst, time_dim='time')

sst_djf  = seasons_sst['DJF']
sst_mam  = seasons_sst['MAM']
sst_jjas = seasons_sst['JJAS']
sst_on   = seasons_sst['ON']


seasons_sss = seasonal_climatology(sss, time_dim='time')

sss_djf  = seasons_sss['DJF']
sss_mam  = seasons_sss['MAM']
sss_jjas = seasons_sss['JJAS']
sss_on   = seasons_sss['ON']

#%% seasonal pco2,sst,sss,dic,alk climatology
# ---- Seasons ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']
xticks = np.arange(35,125,10)
yticks = np.arange(-30,40,10)
# ---- Plot configurations ----
plot_configs_mean = [
    {
        'fields': {
            'DJF': sss_djf,
            'MAM': sss_mam,
            'JJAS': sss_jjas,
            'ON':  sss_on
        },
        'levels': np.arange(30,38,0.25),
        'ticks': np.arange(30,38,0.5),
        'cmap': cmaps.MPL_PuBu,
        'title': 'Sea surface salinity (1994 - 2024)',
        'unit': 'psu'
    },
    {
        'fields': {
            'DJF': sst_djf,
            'MAM': sst_mam,
            'JJAS': sst_jjas,
            'ON':  sst_on
        },
        'levels': np.arange(18,30,0.5),
        'ticks': np.arange(18,30,1),
        'cmap': cmaps.cmp_b2r,
        'title': 'Sea surface temperature (1994 - 2024)',
        'unit': '$^\circ$C'
    },
    {
        'fields': {
            'DJF': tco2_djf,
            'MAM': tco2_mam,
            'JJAS': tco2_jjas,
            'ON':  tco2_on
        },
        'levels': np.arange(1800,2101,10),
        'ticks': np.arange(1800,2101,20),
        'cmap': cmaps.BlGrYeOrReVi200,
        'title': 'Dissolved Inorganic Carbon (1994 - 2024)',
        'unit': 'DIC (micro mol kg-1)'
    },
    {
        'fields': {
            'DJF': talk_djf,
            'MAM': talk_mam,
            'JJAS': talk_jjas,
            'ON':  talk_on
        },
        'levels': np.arange(2100,2401,10),
        'ticks': np.arange(2100,2401,20),
        'cmap': cmaps.BlGrYeOrReVi200,
        'title': 'Alkalinity (1994 - 2024)',
        'unit': 'ALK (micro mol kg-1)'
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

#------- plotting---------

for cfg in plot_configs_mean:

    fig, axs = plt.subplots(
        2, 2, figsize=(18,14),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace':0.06,'hspace': -0.12}
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
        
#-------------contour-------

        # cs = ax.contour(
        #     data.longitude, data.latitude, data,
        #     levels=cfg['levels'],
        #     colors='k',
        #     linewidths=1
        # )
        
        # ax.clabel(cs, inline=True, fontsize= 15, fmt='%d')

        ax.coastlines(zorder =15)
        ax.set_extent([30, 120, -30, 30], crs=ccrs.PlateCarree())
        land = cf.NaturalEarthFeature(
            'physical', 'land', '10m',
            edgecolor='black',
            facecolor='gray'
        )
        ax.add_feature(land, zorder=2)
#-------------title-----
                
        ax.text(0.75, 0.98, f"Mean ({season})",
        transform=ax.transAxes,
        fontsize=18, zorder = 15,fontweight='bold',
        va='top', ha='right')
        # ax.set_title(season, fontsize=18, fontweight='bold')
        
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



        

#-------- Colorbar ----
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks = cfg['ticks'])
    cbar.set_label(cfg['unit'], fontsize= 24,labelpad=13,fontweight = 'bold')
    cbar.ax.tick_params(labelsize=24)

    # plt.suptitle(
    #     f'Seasonal climatology of {cfg["title"]}',
    #     fontsize=28, y=0.87, fontweight = 'bold',
    # )
    
    
    plt.show()

#%% std of pco2,sst,sss,alk,dic 

seasons_pco2 = seasonal_standardeviation(pco2, time_dim='time')

sp_djf_std  = seasons_pco2['DJF']
sp_mam_std  = seasons_pco2['MAM']
sp_jjas_std = seasons_pco2['JJAS']
sp_on_std   = seasons_pco2['ON']
   
seasonal_talk = seasonal_standardeviation(talk, time_dim='time')

talk_djf_std  = seasonal_talk['DJF']
talk_mam_std  = seasonal_talk['MAM']
talk_jjas_std = seasonal_talk['JJAS']
talk_on_std   = seasonal_talk['ON']

seasonal_tco2 = seasonal_standardeviation(tco2, time_dim='time')

tco2_djf_std  = seasonal_tco2['DJF']
tco2_mam_std  = seasonal_tco2['MAM']
tco2_jjas_std = seasonal_tco2['JJAS']
tco2_on_std   = seasonal_tco2['ON']

seasons_sst = seasonal_standardeviation(sst, time_dim='time')

sst_djf_std  = seasons_sst['DJF']
sst_mam_std  = seasons_sst['MAM']
sst_jjas_std = seasons_sst['JJAS']
sst_on_std   = seasons_sst['ON']

seasons_sss = seasonal_standardeviation(sss, time_dim='time')

sss_djf_std  = seasons_sss['DJF']
sss_mam_std  = seasons_sss['MAM']
sss_jjas_std = seasons_sss['JJAS']
sss_on_std   = seasons_sss['ON']

#%% std plot of pco2,sst,sss,dic,alk

# ---- Plot configurations ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']

plot_configs_std = [
    {
        'fields': {
            'DJF': sss_djf_std,
            'MAM': sss_mam_std,
            'JJAS': sss_jjas_std,
            'ON':  sss_on_std
        },
        'levels': 25,
        'cmap': cmaps.MPL_PuBu,
        'title': 'Sea surface salinity (1994 - 2024)',
        'unit': 'PSU'
    },
    {
        'fields': {
            'DJF': sst_djf_std,
            'MAM': sst_mam_std,
            'JJAS': sst_jjas_std,
            'ON':  sst_on_std
        },
        'levels': 25,
        'cmap': cmaps.cmp_b2r,
        'title': 'Sea surface temperature (1994 - 2024)',
        'unit': '$^\circ$C'
    },
    {
        'fields': {
            'DJF': tco2_djf_std,
            'MAM': tco2_mam_std,
            'JJAS': tco2_jjas_std,
            'ON':  tco2_on_std
        },
        'levels': 25,
        'cmap': cmaps.BlGrYeOrReVi200,
        'title': 'Dissolved Inorganic Carbon (1994 - 2024)',
        'unit': 'DIC (micro mol kg-1)'
    },
    {
        'fields': {
            'DJF': talk_djf_std,
            'MAM': talk_mam_std,
            'JJAS': talk_jjas_std,
            'ON':  talk_on_std
        },
        'levels': 25,
        'cmap': cmaps.BlGrYeOrReVi200,
        'title': 'Alkalinity (1994 - 2024)',
        'unit': 'ALK (micro mol kg-1)'
    }
]

# ---- Main plotting loop ----
for cfg in plot_configs_std:

    fig, axs = plt.subplots(
        2, 2, figsize=(18, 14),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace':0.06,'hspace': -0.12}
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
        
#-------------contour

        cs = ax.contour(
            data.longitude, data.latitude, data,
            levels= 10,
            colors='k',
            linewidths=1
        )
        ax.clabel(cs, inline=True, fontsize= 10, fmt='%d')
        
        
        ax.coastlines(zorder =15)
        ax.set_extent([30, 120, -30, 30], crs=ccrs.PlateCarree())
        land = cf.NaturalEarthFeature(
            'physical', 'land', '10m',
            edgecolor='black',
            facecolor='gray'
        )
        ax.add_feature(land, zorder=2)
#-------------title-----
                
        ax.text(0.75, 0.98, f"Mean ({season})",
        transform=ax.transAxes,
        fontsize=18, zorder = 15,fontweight='bold',
        va='top', ha='right')
        # ax.set_title(season, fontsize=18, fontweight='bold')
        
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



        

#-------- Colorbar ----
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(cfg['unit'], fontsize= 24,labelpad=13,fontweight = 'bold')
    cbar.ax.tick_params(labelsize=24)

    # plt.suptitle(
    #     f'Seasonal standard deviation  of {cfg["title"]}',
    #     fontsize=28, y=0.94,fontweight = 'bold',
    # )

    
    plt.show()
    
#%%

for mean_cfg, std_cfg in zip(plot_configs_mean, plot_configs_std):

    fig, axs = plt.subplots(
        2, 4, figsize=(24, 12),
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': 0.08, 'hspace': -0.38}
    )

    # --- reuse land feature (faster & cleaner) ---
    land = cf.NaturalEarthFeature(
        'physical', 'land', '10m',
        edgecolor='black', facecolor='gray'
    )

    # ===================== MEAN =====================
    for i, season in enumerate(seasons):

        ax = axs[i//2, i % 2]
        data = mean_cfg['fields'][season]

        im1 = ax.contourf(
            data.longitude, data.latitude, data,
            levels=mean_cfg['levels'],
            cmap=mean_cfg['cmap'],
            transform=ccrs.PlateCarree(),
            extend='both'
        )

        ax.coastlines(zorder=15)
        ax.set_extent([30, 120, -30, 30])
        ax.set_ylim(-30, 30)
        ax.add_feature(land, zorder=2)

        ax.text(0.75, 0.98, f"Mean ({season})",
                transform=ax.transAxes,
                fontsize=16, fontweight='bold',
                va='top', ha='right')

    # ===================== STD =====================
    for i, season in enumerate(seasons):

        ax = axs[i//2, i % 2 + 2]
        data = std_cfg['fields'][season]

        im2 = ax.contourf(
            data.longitude, data.latitude, data,
            levels=std_cfg['levels'],
            cmap=std_cfg['cmap'],
            transform=ccrs.PlateCarree(),
            extend='both'
        )

        ax.coastlines(zorder=15)
        ax.set_extent([30, 120, -30, 30])
        ax.set_ylim(-30, 30)
        ax.add_feature(land, zorder=2)

        ax.text(0.72, 0.98, f"Std ({season})",
                transform=ax.transAxes,
                fontsize=16, fontweight='bold',
                va='top', ha='right')

    # ===================== TICKS =====================
    for i, ax in enumerate(axs.flat):

        ax.set_xticks(xticks, crs=ccrs.PlateCarree())
        ax.set_yticks(yticks, crs=ccrs.PlateCarree())

        ax.tick_params(
            direction='out',
            length=4,
            width=1,
            labelsize=14
        )

        # hide top row x labels
        if i < 4:
            ax.tick_params(labelbottom=False)

        # hide right side y labels
        if i % 4 != 0:
            ax.tick_params(labelleft=False)

    # ===================== COLORBARS =====================

       # Mean (left block)
    cbar1 = fig.colorbar(
        im1,
        ax=axs[:, :2],
        orientation='horizontal',
        fraction=0.02,   # thinner bar
        pad=0.06,        # closer to plots
        aspect=35        # longer bar
    )
    cbar1.set_label(mean_cfg['unit'], fontsize=16, fontweight='bold')
    cbar1.ax.tick_params(labelsize=13)
    
    # STD (right block)
    cbar2 = fig.colorbar(
        im2,
        ax=axs[:, 2:],
        orientation='horizontal',
        fraction=0.02,
        pad=0.06,
        aspect=35
    )
    cbar2.set_label(std_cfg['unit'], fontsize=16, fontweight='bold')
    cbar2.ax.tick_params(labelsize=13)

    plt.show()
#%% nwio region

nwio_sst = sst.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))

nwio_sss = sss.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))

nwio_pco2 = pco2.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))

nwio_tco2 = tco2.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))

nwio_talk = talk.sel(latitude = slice(5,22.5),longitude = slice(45,65)).mean(dim= ('latitude','longitude'))
#%% nwio_carbonate variables
x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May',
                'Jun','Jul','Aug','Sep','Oct','Nov','Dec']

sst_nwio  = nwio_sst.groupby('time.month').mean('time')
sss_nwio   = nwio_sss.groupby('time.month').mean('time')
pco2_nwio  = nwio_pco2.groupby('time.month').mean('time')
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
ax3.set_ylabel('ALK (µmol kg$^{-1}$)', color='black')
ax3.tick_params(axis='y', labelcolor='black')

# Legend
lines = ln1 + ln2 + ln3
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 0.10),
          ncol=3)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Northwestern Indian Ocean Monthly Climatology of pCO$_2$,DIC,ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%% nwio physical variables
fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# SST (left axis)
ln1 = ax.plot(x, sst_nwio, lw=2, color='navy', label='SST')
ax.set_ylabel('SST (°C)', color='navy')
ax.tick_params(axis='y', labelcolor='navy')

# SSS (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_nwio, lw=2, linestyle='--',
               color='darkgreen', label='SSS')
ax2.set_ylabel('SSS (psu)', color='darkgreen')
ax2.tick_params(axis='y', labelcolor='darkgreen')

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.97, 0.99),
          ncol=2)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Northwestern Indian Ocean Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%% esio region

esio_sst = sst.sel(latitude = slice(-6.6,8),longitude = slice(92,109)).mean(dim= ('latitude','longitude'))

esio_sss = sss.sel(latitude = slice(-6.6,8),longitude = slice(92,109)).mean(dim= ('latitude','longitude'))

esio_pco2 = pco2.sel(latitude = slice(-6.6,8),longitude = slice(92,109)).mean(dim= ('latitude','longitude'))

esio_tco2 = tco2.sel(latitude = slice(-6.6,8),longitude = slice(92,109)).mean(dim= ('latitude','longitude'))

esio_talk = talk.sel(latitude = slice(-6.6,8),longitude = slice(92,109)).mean(dim= ('latitude','longitude'))

#%% esio_carbo 

x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May',
                'Jun','Jul','Aug','Sep','Oct','Nov','Dec']

sst_esio  = esio_sst.groupby('time.month').mean('time')
sss_esio  = esio_sss.groupby('time.month').mean('time')
pco2_esio = esio_pco2.groupby('time.month').mean('time')
dic_esio  = esio_tco2.groupby('time.month').mean('time')
alk_esio  = esio_talk.groupby('time.month').mean('time')


fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# pCO2 (left axis)
ln1 = ax.plot(x, pco2_esio, lw=2, color='firebrick', label='pCO$_2$')
ax.set_ylabel('pCO$_2$ (µatm)', color='firebrick')
ax.tick_params(axis='y', labelcolor='firebrick')

# DIC (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, dic_esio, lw=2, color='purple', label='DIC',linestyle='--')
ax2.set_ylabel('DIC (µmol kg$^{-1}$)', color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# ALK (third axis)
ax3 = ax.twinx()
ax3.spines['right'].set_position(('outward', 60))
ln3 = ax3.plot(x, alk_esio, lw=2, color='black', label='ALK',linestyle='--')
ax3.set_ylabel('ALK (µmol kg$^{-1}$)', color='black')
ax3.tick_params(axis='y', labelcolor='black')

# Legend
lines = ln1 + ln2 + ln3
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 0.99),
          ncol=3)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Eastern Indian Ocean Monthly Climatology of pCO$_2$,DIC,ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%% esio_phys
fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# SST (left axis)
ln1 = ax.plot(x, sst_esio, lw=2, color='navy', label='SST')
ax.set_ylabel('SST (°C)', color='navy')
ax.tick_params(axis='y', labelcolor='navy')

# SSS (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_esio, lw=2, linestyle='--',
               color='darkgreen', label='SSS')
ax2.set_ylabel('SSS (psu)', color='darkgreen')
ax2.tick_params(axis='y', labelcolor='darkgreen')

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 0.99),
          ncol=2)
# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Eastern Indian Ocean Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()

#%% eio region

eio_sst = sst.sel(latitude = slice(-6.5,5),longitude = slice(49,92)).mean(dim= ('latitude','longitude'))

eio_sss = sss.sel(latitude = slice(-6.5,5),longitude = slice(49,92)).mean(dim= ('latitude','longitude'))

eio_pco2 = pco2.sel(latitude = slice(-6.5,5),longitude = slice(49,92)).mean(dim= ('latitude','longitude'))

eio_tco2 = tco2.sel(latitude = slice(-6.5,5),longitude = slice(49,92)).mean(dim= ('latitude','longitude'))

eio_talk = talk.sel(latitude = slice(-6.5,5),longitude = slice(49,92)).mean(dim= ('latitude','longitude'))

#%% eio_carbo 

x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May',
                'Jun','Jul','Aug','Sep','Oct','Nov','Dec']

sst_eio  = eio_sst.groupby('time.month').mean('time')
sss_eio  = eio_sss.groupby('time.month').mean('time')
pco2_eio = eio_pco2.groupby('time.month').mean('time')
dic_eio  = eio_tco2.groupby('time.month').mean('time')
alk_eio  = eio_talk.groupby('time.month').mean('time')


fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# pCO2 (left axis)
ln1 = ax.plot(x, pco2_eio, lw=2, color='firebrick', label='pCO$_2$')
ax.set_ylabel('pCO$_2$ (µatm)', color='firebrick')
ax.tick_params(axis='y', labelcolor='firebrick')

# DIC (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, dic_eio, lw=2, color='purple', label='DIC',linestyle='--')
ax2.set_ylabel('DIC (µmol kg$^{-1}$)', color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# ALK (third axis)
ax3 = ax.twinx()
ax3.spines['right'].set_position(('outward', 60))
ln3 = ax3.plot(x, alk_eio, lw=2, color='black', label='ALK',linestyle='--')
ax3.set_ylabel('ALK (µmol kg$^{-1}$)', color='black')
ax3.tick_params(axis='y', labelcolor='black')

# Legend
lines = ln1 + ln2 + ln3
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 0.10),
          ncol=3)
# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Equatorial Indian Ocean Monthly Climatology of pCO$_2$,DIC,ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
# eio_phys
fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# SST (left axis)
ln1 = ax.plot(x, sst_eio, lw=2, color='navy', label='SST')
ax.set_ylabel('SST (°C)', color='navy')
ax.tick_params(axis='y', labelcolor='navy')

# SSS (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_eio, lw=2, linestyle='--',
               color='darkgreen', label='SSS')
ax2.set_ylabel('SSS (psu)', color='darkgreen')
ax2.tick_params(axis='y', labelcolor='darkgreen')

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 0.10),
          ncol=2)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Equatorial Indian Ocean Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()

#%% nas

nas_sst = sst.sel(latitude = slice(22.5, 28),longitude = slice(56, 70)).mean(dim= ('latitude','longitude'))

nas_sss = sss.sel(latitude = slice(22.5, 28),longitude = slice(56, 70)).mean(dim= ('latitude','longitude'))

nas_pco2 = pco2.sel(latitude = slice(22.5, 28),longitude = slice(56, 70)).mean(dim= ('latitude','longitude'))

nas_tco2 = tco2.sel(latitude = slice(22.5, 28),longitude = slice(56, 70)).mean(dim= ('latitude','longitude'))

nas_talk = talk.sel(latitude = slice(22.5, 28),longitude = slice(56, 70)).mean(dim= ('latitude','longitude'))

#%% nas_carbo

x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May',
                'Jun','Jul','Aug','Sep','Oct','Nov','Dec']

sst_nas  = nas_sst.groupby('time.month').mean('time')
sss_nas  = nas_sss.groupby('time.month').mean('time')
pco2_nas = nas_pco2.groupby('time.month').mean('time')
dic_nas  = nas_tco2.groupby('time.month').mean('time')
alk_nas  = nas_talk.groupby('time.month').mean('time')


fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# pCO2 (left axis)
ln1 = ax.plot(x, pco2_nas, lw=2, color='firebrick', label='pCO$_2$')
ax.set_ylabel('pCO$_2$ (µatm)', color='firebrick')
ax.tick_params(axis='y', labelcolor='firebrick')

# DIC (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, dic_nas, lw=2, color='purple', label='DIC',linestyle='--')
ax2.set_ylabel('DIC (µmol kg$^{-1}$)', color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# ALK (third axis)
ax3 = ax.twinx()
ax3.spines['right'].set_position(('outward', 60))
ln3 = ax3.plot(x, alk_nas, lw=2, color='black', label='ALK',linestyle='--')
ax3.set_ylabel('ALK (µmol kg$^{-1}$)', color='black')
ax3.tick_params(axis='y', labelcolor='black')

# Legend
lines = ln1 + ln2 + ln3
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 0.10),
          ncol=3)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Northern Arabian Sea Monthly Climatology of pCO$_2$,DIC,ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%% nas_phys
fig, ax = plt.subplots(figsize=(11, 5),dpi = 600)

# SST (left axis)
ln1 = ax.plot(x, sst_nas, lw=2, color='navy', label='SST')
ax.set_ylabel('SST (°C)', color='navy')
ax.tick_params(axis='y', labelcolor='navy')

# SSS (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_nas, lw=2, linestyle='--',
               color='darkgreen', label='SSS')
ax2.set_ylabel('SSS (psu)', color='darkgreen')
ax2.tick_params(axis='y', labelcolor='darkgreen')

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 0.10),
          ncol=2)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Northern Arabian Sea Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%%


pco2_clr= '#F40009'
dic_clr = 'darkred'
alk_clr = 'b'
sst_clr = '#BE8500'
sss_clr = '#0039BE'

#%% nwio carbonate(edited)
x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May',
                'Jun','Jul','Aug','Sep','Oct','Nov','Dec']


t_mean    = nwio_sst.groupby('time.month').mean('time')
t_std     = nwio_sst.groupby('time.month').std('time')

s_mean    = nwio_sss.groupby('time.month').mean('time')
s_std     = nwio_sss.groupby('time.month').std('time')

pco2_mean = nwio_pco2.groupby('time.month').mean('time')
pco2_std  = nwio_pco2.groupby('time.month').std('time')

dic_mean  = nwio_tco2.groupby('time.month').mean('time')
dic_std   = nwio_tco2.groupby('time.month').std('time')

alk_mean  = nwio_talk.groupby('time.month').mean('time')
alk_std   = nwio_talk.groupby('time.month').std('time')

fig, ax = plt.subplots(figsize=(13,5),dpi = 600)

# DIC (left axis)
ln1 = ax.plot(x, dic_mean, lw=2, color=dic_clr, label='DIC')

ax.fill_between(
    x,
    dic_mean - dic_std,
    dic_mean + dic_std,
    color=dic_clr,
    alpha=0.25,
    linewidth=3
)

ax.set_ylabel('DIC (µmol kg$^{-1}$)', color=dic_clr)
ax.tick_params(axis='y', labelcolor=dic_clr)

# ALK (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, alk_mean, lw=2, color=alk_clr,
               linestyle='--', label='ALK')

ax2.fill_between(
    x,
    alk_mean - alk_std,
    alk_mean + alk_std,
    color=alk_clr,
    alpha=0.2,
    linewidth=3
)

ax2.set_ylabel('ALK (µmol kg$^{-1}$)', color=alk_clr)
ax2.tick_params(axis='y', labelcolor=alk_clr)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 0.10),
          ncol=2)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Northwestern Indian Ocean Monthly Climatology of DIC & ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()


#%% nwio_phy (edited)

fig, ax = plt.subplots(figsize=(13, 5),dpi = 600)

# SST (left axis)
ln1 = ax.plot(x, t_mean, lw=2, color=sst_clr, label='SST')

ax.fill_between(
    x,
    t_mean - t_std,
    t_mean + t_std,
    color=sst_clr,
    alpha=0.2,linewidth=3
)
ax.set_ylabel('SST (°C)', color=sst_clr)
ax.tick_params(axis='y', labelcolor=sst_clr)

# SSS (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, s_mean, lw=2, linestyle='--',
               color=sss_clr, label='SSS')

ax2.fill_between(
    x,
    s_mean - s_std,
    s_mean + s_std,
    color=sss_clr,
    alpha=0.2,linewidth=3
)

ax2.set_ylabel('SSS (psu)', color=sss_clr)
ax2.tick_params(axis='y', labelcolor=sss_clr)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='upper right',
    bbox_to_anchor=(0.99,0.10),
    ncol=3)
# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Northwestern Indian Ocean Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%% esio_carbo_edited
# ==============================
# Monthly climatology (MEAN)
# ==============================
sst_mean  = esio_sst.groupby('time.month').mean('time')
sss_mean  = esio_sss.groupby('time.month').mean('time')
pco2_mean = esio_pco2.groupby('time.month').mean('time')
dic_mean  = esio_tco2.groupby('time.month').mean('time')
alk_mean  = esio_talk.groupby('time.month').mean('time')

# ==============================
# Monthly STD (interannual variability)
# ==============================
sst_std  = esio_sst.groupby('time.month').std('time')
sss_std  = esio_sss.groupby('time.month').std('time')
pco2_std = esio_pco2.groupby('time.month').std('time')
dic_std  = esio_tco2.groupby('time.month').std('time')
alk_std  = esio_talk.groupby('time.month').std('time')

x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# ==============================
# Carbon system plot (with STD)
# ==============================
fig, ax = plt.subplots(figsize=(13,5),dpi = 600)

# DIC (left axis)
ln1 = ax.plot(x, dic_mean, lw=2, color=dic_clr, label='DIC')

ax.fill_between(
    x,
    dic_mean - dic_std,
    dic_mean + dic_std,
    color=dic_clr,
    alpha=0.2,
    linewidth=3
)

ax.set_ylabel('DIC (µmol kg$^{-1}$)', color=dic_clr)
ax.tick_params(axis='y', labelcolor=dic_clr)

# ALK (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, alk_mean, lw=2, color=alk_clr,
               linestyle='--', label='ALK')

ax2.fill_between(
    x,
    alk_mean - alk_std,
    alk_mean + alk_std,
    color=alk_clr,
    alpha=0.2,
    linewidth=3
)

ax2.set_ylabel('ALK (µmol kg$^{-1}$)', color=alk_clr)
ax2.tick_params(axis='y', labelcolor=alk_clr)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='upper right', ncol=2)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Eastern Indian Ocean Monthly Climatology of DIC & ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()


#%% esio_phy_edited

fig, ax = plt.subplots(figsize=(13,5), dpi=600)

# SST
ln1 = ax.plot(x, sst_mean, lw=2, color=sst_clr, label='SST')
ax.fill_between(x, sst_mean - sst_std, sst_mean + sst_std,
                color=sst_clr, alpha=0.2,linewidth=3)

ax.set_ylabel('SST (°C)', color=sst_clr)
ax.tick_params(axis='y', labelcolor=sst_clr)

# SSS
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_mean, lw=2, linestyle='--',
               color=sss_clr, label='SSS')
ax2.fill_between(x, sss_mean - sss_std, sss_mean + sss_std,
                 color=sss_clr, alpha=0.2,linewidth=3)

ax2.set_ylabel('SSS (psu)', color=sss_clr)
ax2.tick_params(axis='y', labelcolor=sss_clr)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,loc='upper right',ncol =3)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Eastern Indian Ocean Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%% eio_carbo_edited

# ==============================
# Monthly climatology (MEAN)
# ==============================
sst_mean  = eio_sst.groupby('time.month').mean('time')
sss_mean  = eio_sss.groupby('time.month').mean('time')
pco2_mean = eio_pco2.groupby('time.month').mean('time')
dic_mean  = eio_tco2.groupby('time.month').mean('time')
alk_mean  = eio_talk.groupby('time.month').mean('time')

# ==============================
# Monthly STD (interannual variability)
# ==============================
sst_std  = eio_sst.groupby('time.month').std('time')
sss_std  = eio_sss.groupby('time.month').std('time')
pco2_std = eio_pco2.groupby('time.month').std('time')
dic_std  = eio_tco2.groupby('time.month').std('time')
alk_std  = eio_talk.groupby('time.month').std('time')

# Convert to numpy (safe plotting)
sst_mean, sst_std = sst_mean.values, sst_std.values
sss_mean, sss_std = sss_mean.values, sss_std.values
pco2_mean, pco2_std = pco2_mean.values, pco2_std.values
dic_mean, dic_std = dic_mean.values, dic_std.values
alk_mean, alk_std = alk_mean.values, alk_std.values

# ==============================
# Common axis
# ==============================
x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# ==============================
# Carbon system plot
# ==============================
fig, ax = plt.subplots(figsize=(13,5), dpi=600)

# DIC (left axis)
ln1 = ax.plot(x, dic_mean, lw=2, color=dic_clr, label='DIC')

ax.fill_between(
    x,
    dic_mean - dic_std,
    dic_mean + dic_std,
    color=dic_clr,
    alpha=0.2,
    linewidth=3
)

ax.set_ylabel('DIC (µmol kg$^{-1}$)', color=dic_clr)
ax.tick_params(axis='y', labelcolor=dic_clr)

# ALK (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, alk_mean, lw=2, color=alk_clr,
               linestyle='--', label='ALK')

ax2.fill_between(
    x,
    alk_mean - alk_std,
    alk_mean + alk_std,
    color=alk_clr,
    alpha=0.2,
    linewidth=3
)

ax2.set_ylabel('ALK (µmol kg$^{-1}$)', color=alk_clr)
ax2.tick_params(axis='y', labelcolor=alk_clr)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='best',
          bbox_to_anchor=(0.99, 0.10),
          ncol=2)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Equatorial Indian Ocean Monthly Climatology of DIC & ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()



# ==============================
# Physical variables plot (SST, SSS)
# ==============================
fig, ax = plt.subplots(figsize=(13,5), dpi=600)

# SST
ln1 = ax.plot(x, sst_mean, lw=2, color=sst_clr, label='SST')
ax.fill_between(x, sst_mean - sst_std, sst_mean + sst_std,
                color=sst_clr, alpha=0.2,linewidth=3)

ax.set_ylabel('SST (°C)', color=sst_clr)
ax.tick_params(axis='y', labelcolor=sst_clr)

# SSS
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_mean, lw=2, linestyle='--',
               color=sss_clr, label='SSS')

ax2.fill_between(x, sss_mean - sss_std, sss_mean + sss_std,
                 color=sss_clr, alpha=0.2,linewidth=3)

ax2.set_ylabel('SSS (psu)', color=sss_clr)
ax2.tick_params(axis='y', labelcolor=sss_clr)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='upper left',ncol =3)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Equatorial Indian Ocean Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%% nas_carbo_phy_edited
# ==============================
# Monthly climatology (MEAN)
# ==============================
sst_mean  = nas_sst.groupby('time.month').mean('time')
sss_mean  = nas_sss.groupby('time.month').mean('time')
pco2_mean = nas_pco2.groupby('time.month').mean('time')
dic_mean  = nas_tco2.groupby('time.month').mean('time')
alk_mean  = nas_talk.groupby('time.month').mean('time')

# ==============================
# Monthly STD (interannual variability)
# ==============================
sst_std  = nas_sst.groupby('time.month').std('time')
sss_std  = nas_sss.groupby('time.month').std('time')
pco2_std = nas_pco2.groupby('time.month').std('time')
dic_std  = nas_tco2.groupby('time.month').std('time')
alk_std  = nas_talk.groupby('time.month').std('time')

# Convert to numpy (safe plotting)
sst_mean, sst_std = sst_mean.values, sst_std.values
sss_mean, sss_std = sss_mean.values, sss_std.values
pco2_mean, pco2_std = pco2_mean.values, pco2_std.values
dic_mean, dic_std = dic_mean.values, dic_std.values
alk_mean, alk_std = alk_mean.values, alk_std.values

# ==============================
# Common axis
# ==============================
x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# ==============================
# Carbon system plot
# ==============================
fig, ax = plt.subplots(figsize=(13,5), dpi=600)

# DIC (left axis)
ln1 = ax.plot(x, dic_mean, lw=2, color=dic_clr, label='DIC')

ax.fill_between(
    x,
    dic_mean - dic_std,
    dic_mean + dic_std,
    color=dic_clr,
    alpha=0.2,
    linewidth=3
)

ax.set_ylabel('DIC (µmol kg$^{-1}$)', color=dic_clr)
ax.tick_params(axis='y', labelcolor=dic_clr)

# ALK (right axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, alk_mean, lw=2, color=alk_clr,
               linestyle='--', label='ALK')

ax2.fill_between(
    x,
    alk_mean - alk_std,
    alk_mean + alk_std,
    color=alk_clr,
    alpha=0.2,
    linewidth=3
)

ax2.set_ylabel('ALK (µmol kg$^{-1}$)', color=alk_clr)
ax2.tick_params(axis='y', labelcolor=alk_clr)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='best',
          bbox_to_anchor=(0.99, 0.10),
          ncol=2)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Northern Arabian Sea Monthly Climatology of DIC & ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
# ==============================
# Physical variables plot (SST, SSS)
# ==============================
fig, ax = plt.subplots(figsize=(13,5), dpi=600)

# SST
ln1 = ax.plot(x, sst_mean, lw=2, color=sst_clr, label='SST')
ax.fill_between(x, sst_mean - sst_std, sst_mean + sst_std,
                color=sst_clr, alpha=0.2,linewidth=3)

ax.set_ylabel('SST (°C)', color=sst_clr)
ax.tick_params(axis='y', labelcolor=sst_clr)

# SSS
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_mean, lw=2, linestyle='--',
               color=sss_clr, label='SSS')

ax2.fill_between(x, sss_mean - sss_std, sss_mean + sss_std,
                 color=sss_clr, alpha=0.2,linewidth=3)

ax2.set_ylabel('SSS (psu)', color=sss_clr)
ax2.tick_params(axis='y', labelcolor=sss_clr)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='best',ncol = 3)

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels)
ax.set_xlim(0, 11)
ax.set_xlabel('Month')
# ax.set_title('Northern Arabian Sea Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()

#%% pco2 _subplot

# ==============================
# Monthly climatology (MEAN & STD)
# ==============================

regions = {
    'NWIO': nwio_pco2,
    'NAS' : nas_pco2,
    'ESIO': esio_pco2,
    'EIO' : eio_pco2
}

# Common axis
x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May','Jun',
                'Jul','Aug','Sep','Oct','Nov','Dec']

# ==============================
# Plot
# ==============================
fig, axs = plt.subplots(2, 2, figsize=(14,10), dpi=600)

for i, (name, data) in enumerate(regions.items()):

    ax = axs.flat[i]

    # Monthly mean & std
    mean = data.groupby('time.month').mean('time')
    std  = data.groupby('time.month').std('time')

    # Convert to numpy (safe plotting)
    mean = mean.values
    std  = std.values

    # Plot
    ax.plot(x, mean, lw=2, color=pco2_clr, label='pCO$_2$')

    ax.fill_between(
        x,
        mean - std,
        mean + std,
        color=pco2_clr,
        alpha=0.25,
        linewidth=2
    )

    # Formatting
    ax.set_title(name, fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(month_labels)
    ax.set_xlim(0, 11)
    ax.set_ylim(340,460)  # adjust based on your data
    ax.grid(True, linestyle='--', alpha=0.3)

    # Y label only on left column
    if i % 2 == 0:
        ax.set_ylabel('pCO$_2$ (µatm)', fontsize=14)

    # X label only on bottom row
    if i >= 2:
        ax.set_xlabel('Month', fontsize=14)

# ==============================
# Super title
# ==============================
# plt.suptitle(
#     'Monthly Climatology of pCO$_2$',
#     fontsize=20, fontweight='bold', y=0.95
# )

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.show()