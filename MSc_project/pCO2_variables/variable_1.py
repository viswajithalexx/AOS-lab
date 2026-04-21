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


file1 = '/home/bobco-08/24cl05012/CO2/data/data_1/cmems/bio_var/cmems_obs-mob_glo_bgc-car_my_irr-i_1774948756844_30N-30S_30E-120E_1985-2024.nc'

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
#     {
#         'fields': {
#             'DJF': sss_djf,
#             'MAM': sss_mam,
#             'JJAS': sss_jjas,
#             'ON':  sss_on
#         },
#         'levels': np.arange(30,38,0.25),
#         'ticks': np.arange(30,38,0.5),
#         'cmap': cmaps.MPL_PuBu,
#         'title': 'Sea surface salinity (1994 - 2024)',
#         'unit': 'psu'
#     },
    {
        'fields': {
            'DJF': sst_djf,
            'MAM': sst_mam,
            'JJAS': sst_jjas,
            'ON':  sst_on
        },
        'levels': np.arange(22,30.2,0.2),
        'ticks': np.arange(22,31,1),
        'cmap': cmaps.cmp_b2r,
        'title': 'Sea surface temperature (1994 - 2024)',
        'unit': '$^\circ$C'
    },
    # {
    #     'fields': {
    #         'DJF': tco2_djf,
    #         'MAM': tco2_mam,
    #         'JJAS': tco2_jjas,
    #         'ON':  tco2_on
    #     },
    #     'levels': np.arange(1800,2101,10),
    #     'ticks': np.arange(1800,2101,20),
    #     'cmap': cmaps.BlGrYeOrReVi200,
    #     'title': 'Dissolved Inorganic Carbon (1994 - 2024)',
    #     'unit': 'DIC (micro mol kg-1)'
    # },
    # {
    #     'fields': {
    #         'DJF': talk_djf,
    #         'MAM': talk_mam,
    #         'JJAS': talk_jjas,
    #         'ON':  talk_on
    #     },
    #     'levels': np.arange(2100,2401,10),
    #     'ticks': np.arange(2100,2401,20),
    #     'cmap': cmaps.BlGrYeOrReVi200,
    #     'title': 'Alkalinity (1994 - 2024)',
    #     'unit': 'ALK (micro mol kg-1)'
    # }
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
        2, 2, figsize=(18,14),dpi = 600,
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
        #     levels=25,
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

months = np.arange(1, 13)

print("\n===== Monthly Climatology (NWIO) =====\n")

for m in months:
    print(f"Month {m:02d}:")
    print(f"  SST  : {sst_nwio.sel(month=m).values:.2f}")
    print(f"  SSS  : {sss_nwio.sel(month=m).values:.2f}")
    print(f"  pCO2 : {pco2_nwio.sel(month=m).values:.2f}")
    print(f"  DIC  : {dic_nwio.sel(month=m).values:.2f}")
    print(f"  ALK  : {alk_nwio.sel(month=m).values:.2f}")
    print("-" * 40)
#%%
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

fig, ax = plt.subplots(figsize=(13,5))

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
ax.set_title('Northwestern Indian Ocean Monthly Climatology of DIC & ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()


#%% nwio_phy (edited)

fig, ax = plt.subplots(figsize=(13, 5))

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
ax.set_title('Northwestern Indian Ocean Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()

#%% subplot 1 by 2 nwio
fig, axs = plt.subplots(2, 1, figsize=(10, 9),dpi = 300,
                        gridspec_kw={'hspace': 0.3})

# Global font sizes
label_fs = 14
tick_fs = 13
title_fs = 14
legend_fs = 9

# =========================
# LEFT PANEL: DIC & ALK
# =========================
ax = axs[0]

# DIC
ln1 = ax.plot(x, dic_mean, lw=2, color=dic_clr, label='DIC')
ax.fill_between(x, dic_mean - dic_std, dic_mean + dic_std,
                color=dic_clr, alpha=0.25)

ax.set_ylabel('DIC (µmol kg$^{-1}$)', color=dic_clr, fontsize=label_fs)
ax.tick_params(axis='y', labelcolor=dic_clr, labelsize=tick_fs)

# ALK (twin axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, alk_mean, lw=2, linestyle='--',
               color=alk_clr, label='ALK')
ax2.fill_between(x, alk_mean - alk_std, alk_mean + alk_std,
                 color=alk_clr, alpha=0.2)

ax2.set_ylabel('ALK (µmol kg$^{-1}$)', color=alk_clr, fontsize=label_fs)
ax2.tick_params(axis='y', labelcolor=alk_clr, labelsize=tick_fs)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.35, 1.0),  # (x, y)
          fontsize=legend_fs,
          ncol=len(labels))

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=tick_fs)
ax.set_xlim(0, 11)
ax.set_title('NWIO: DIC & ALK', fontsize=title_fs, fontweight ='bold')
ax.grid(True, linestyle='--', alpha=0.3)

# =========================
# RIGHT PANEL: SST & SSS
# =========================
ax = axs[1]

# SST
ln1 = ax.plot(x, t_mean, lw=2, color=sst_clr, label='SST')
ax.fill_between(x, t_mean - t_std, t_mean + t_std,
                color=sst_clr, alpha=0.2)

ax.set_ylabel('SST (°C)', color=sst_clr, fontsize=label_fs)
ax.tick_params(axis='y', labelcolor=sst_clr, labelsize=tick_fs)

# SSS (twin axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, s_mean, lw=2, linestyle='--',
               color=sss_clr, label='SSS')
ax2.fill_between(x, s_mean - s_std, s_mean + s_std,
                 color=sss_clr, alpha=0.2)

ax2.set_ylabel('SSS (psu)', color=sss_clr, fontsize=label_fs)
ax2.tick_params(axis='y', labelcolor=sss_clr, labelsize=tick_fs)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.97, 1.0),  # (x, y)
          fontsize=legend_fs,
          ncol=len(labels))

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=tick_fs)
ax.set_xlim(0, 11)
ax.set_xlabel('Month', fontsize=label_fs)
ax.set_title('NWIO: SST & SSS', fontsize=title_fs,fontweight ='bold')
ax.grid(True, linestyle='--', alpha=0.3)

# =========================
# FINAL LAYOUT
# =========================
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
fig, ax = plt.subplots(figsize=(13,5))

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
ax.set_title('Eastern Indian Ocean Monthly Climatology of DIC & ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()


#%% esio_phy_edited

fig, ax = plt.subplots(figsize=(13,5))

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
ax.set_title('Eastern Indian Ocean Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()

#%% subplot 1 by 2 esio

fig, axs = plt.subplots(2, 1, figsize=(10, 9),dpi = 300,
                        gridspec_kw={'hspace': 0.3})
#Global font sizes
label_fs = 15
tick_fs = 13
title_fs = 15
legend_fs = 9

# =========================
# LEFT PANEL: DIC & ALK
# =========================
ax = axs[0]

# DIC
ln1 = ax.plot(x, dic_mean, lw=2, color=dic_clr, label='DIC')
ax.fill_between(x, dic_mean - dic_std, dic_mean + dic_std,
                color=dic_clr, alpha=0.25)

ax.set_ylabel('DIC (µmol kg$^{-1}$)', color=dic_clr, fontsize=label_fs)
ax.tick_params(axis='y', labelcolor=dic_clr, labelsize=tick_fs)

# ALK (twin axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, alk_mean, lw=2, linestyle='--',
               color=alk_clr, label='ALK')
ax2.fill_between(x, alk_mean - alk_std, alk_mean + alk_std,
                 color=alk_clr, alpha=0.2)

ax2.set_ylabel('ALK (µmol kg$^{-1}$)', color=alk_clr, fontsize=label_fs)
ax2.tick_params(axis='y', labelcolor=alk_clr, labelsize=tick_fs)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 1.0),  # (x, y)
          fontsize=legend_fs,
          ncol=len(labels))

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=tick_fs)
ax.set_xlim(0, 11)
ax.set_title('ESIO: DIC & ALK', fontsize=title_fs, fontweight ='bold')
ax.grid(True, linestyle='--', alpha=0.3)

# =========================
# RIGHT PANEL: SST & SSS
# =========================
ax = axs[1]

# SST
ln1 = ax.plot(x, sst_mean, lw=2, color=sst_clr, label='SST')
ax.fill_between(x, sst_mean - t_std, sst_mean + t_std,
                color=sst_clr, alpha=0.2)

ax.set_ylabel('SST (°C)', color=sst_clr, fontsize=label_fs)
ax.tick_params(axis='y', labelcolor=sst_clr, labelsize=tick_fs)

# SSS (twin axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_mean, lw=2, linestyle='--',
               color=sss_clr, label='SSS')
ax2.fill_between(x, sss_mean - s_std, sss_mean + s_std,
                 color=sss_clr, alpha=0.2)

ax2.set_ylabel('SSS (psu)', color=sss_clr, fontsize=label_fs)
ax2.tick_params(axis='y', labelcolor=sss_clr, labelsize=tick_fs)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.99, 1.0),  # (x, y)
          fontsize=legend_fs,
          ncol=len(labels))

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=tick_fs)
ax.set_xlim(0, 11)
ax.set_xlabel('Month', fontsize=label_fs)
ax.set_title('ESIO: SST & SSS', fontsize=title_fs,fontweight ='bold')
ax.grid(True, linestyle='--', alpha=0.3)

# =========================
# FINAL LAYOUT
# =========================
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
fig, ax = plt.subplots(figsize=(13,5))

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
ax.set_title('Equatorial Indian Ocean Monthly Climatology of DIC & ALK')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()



# ==============================
# Physical variables plot (SST, SSS)
# ==============================
fig, ax = plt.subplots(figsize=(13,5))

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
ax.set_title('Equatorial Indian Ocean Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()

#%% subplot 1 by 2

fig, axs = plt.subplots(2, 1, figsize=(10, 9),dpi =300,
                        gridspec_kw={'hspace': 0.3})

# =========================
# LEFT PANEL: DIC & ALK
# =========================
ax = axs[0]

# DIC
ln1 = ax.plot(x, dic_mean, lw=2, color=dic_clr, label='DIC')
ax.fill_between(x, dic_mean - dic_std, dic_mean + dic_std,
                color=dic_clr, alpha=0.25)


ax.set_ylabel('DIC (µmol kg$^{-1}$)', color=dic_clr, fontsize=label_fs)
ax.tick_params(axis='y', labelcolor=dic_clr, labelsize=tick_fs)


# ALK (twin axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, alk_mean, lw=2, linestyle='--',
               color=alk_clr, label='ALK')
ax2.fill_between(x, alk_mean - alk_std, alk_mean + alk_std,
                 color=alk_clr, alpha=0.2)

ax2.set_ylabel('ALK (µmol kg$^{-1}$)', color=alk_clr, fontsize=label_fs)
ax2.tick_params(axis='y', labelcolor=alk_clr, labelsize=tick_fs)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.21, 1.0),  # (x, y)
          fontsize=legend_fs,
          ncol=len(labels))

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=tick_fs)
ax.set_xlim(0, 11)
ax.set_title('EIO: DIC & ALK', fontsize=title_fs, fontweight ='bold')
ax.grid(True, linestyle='--', alpha=0.3)
# =========================
# RIGHT PANEL: SST & SSS
# =========================
ax = axs[1]

# SST
ln1 = ax.plot(x, sst_mean, lw=2, color=sst_clr, label='SST')
ax.fill_between(x, sst_mean - t_std, sst_mean + t_std,
                color=sst_clr, alpha=0.2)

ax.set_ylabel('SST (°C)', color=sst_clr, fontsize=label_fs)
ax.tick_params(axis='y', labelcolor=sst_clr, labelsize=tick_fs)

# SSS (twin axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_mean, lw=2, linestyle='--',
               color=sss_clr, label='SSS')
ax2.fill_between(x, sss_mean - s_std, sss_mean + s_std,
                 color=sss_clr, alpha=0.2)

ax2.set_ylabel('SSS (psu)', color=sss_clr, fontsize=label_fs)
ax2.tick_params(axis='y', labelcolor=sss_clr, labelsize=tick_fs)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.21, 1.0),  # (x, y)
          fontsize=legend_fs,
          ncol=len(labels))

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=tick_fs)
ax.set_xlim(0, 11)
ax.set_xlabel('Month', fontsize=label_fs)
ax.set_title('EIO: SST & SSS', fontsize=title_fs,fontweight ='bold')
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
fig, ax = plt.subplots(figsize=(13,5))

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
fig, ax = plt.subplots(figsize=(13,5))

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
ax.set_title('Northern Arabian Sea Monthly Climatology of SST,SSS')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
#%% subplot 1 by 2
fig, axs = plt.subplots(2, 1, figsize=(10, 9),dpi = 300,
                        gridspec_kw={'hspace': 0.3})

# =========================
# LEFT PANEL: DIC & ALK
# =========================
ax = axs[0]

# DIC
ln1 = ax.plot(x, dic_mean, lw=2, color=dic_clr, label='DIC')
ax.fill_between(x, dic_mean - dic_std, dic_mean + dic_std,
                color=dic_clr, alpha=0.25)

ax.set_ylabel('DIC (µmol kg$^{-1}$)', color=dic_clr, fontsize=label_fs)
ax.tick_params(axis='y', labelcolor=dic_clr, labelsize=tick_fs)

# ALK (twin axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, alk_mean, lw=2, linestyle='--',
               color=alk_clr, label='ALK')
ax2.fill_between(x, alk_mean - alk_std, alk_mean + alk_std,
                 color=alk_clr, alpha=0.2)

ax2.set_ylabel('ALK (µmol kg$^{-1}$)', color=alk_clr, fontsize=label_fs)
ax2.tick_params(axis='y', labelcolor=alk_clr, labelsize=tick_fs)


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
ax.set_title('NAS: DIC & ALK',fontsize=title_fs, fontweight ='bold')
ax.grid(True, linestyle='--', alpha=0.3)

# =========================
# RIGHT PANEL: SST & SSS
# =========================
ax = axs[1]

# SST
ln1 = ax.plot(x, sst_mean, lw=2, color=sst_clr, label='SST')
ax.fill_between(x, sst_mean - t_std, sst_mean + t_std,
                color=sst_clr, alpha=0.2)

ax.set_ylabel('SST (°C)', color=sst_clr, fontsize=label_fs)
ax.tick_params(axis='y', labelcolor=sst_clr, labelsize=tick_fs)

# SSS (twin axis)
ax2 = ax.twinx()
ln2 = ax2.plot(x, sss_mean, lw=2, linestyle='--',
               color=sss_clr, label='SSS')
ax2.fill_between(x, sss_mean - s_std, sss_mean + s_std,
                 color=sss_clr, alpha=0.2)

ax2.set_ylabel('SSS (psu)', color=sss_clr, fontsize=label_fs)
ax2.tick_params(axis='y', labelcolor=sss_clr, labelsize=tick_fs)

# Legend
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,
          loc='upper right',
          bbox_to_anchor=(0.22, 1.0),  # (x, y)
          fontsize=legend_fs,
          ncol=len(labels))

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(month_labels, fontsize=tick_fs)
ax.set_xlim(0, 11)
ax.set_xlabel('Month', fontsize=label_fs)
ax.set_title('NAS: SST & SSS', fontsize=title_fs,fontweight ='bold')
ax.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.show()
# =========================
# FINAL LAYOUT
# =========================
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
fig, axs = plt.subplots(2, 2, figsize=(14,10))

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
    # ax.set_ylim(340,460)  # adjust based on your data
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
#%%
import numpy as np
import matplotlib.pyplot as plt

regions = [
    ('NWIO', nwio_sst, nwio_sss, nwio_tco2, nwio_talk),
    ('ESIO', esio_sst, esio_sss, esio_tco2, esio_talk),
    ('EIO',  eio_sst,  eio_sss,  eio_tco2,  eio_talk),
    ('NAS',  nas_sst,  nas_sss,  nas_tco2,  nas_talk),
]

x = np.arange(12)
month_labels = ['Jan','Feb','Mar','Apr','May','Jun',
                'Jul','Aug','Sep','Oct','Nov','Dec']

fig, axes = plt.subplots(len(regions), 2, figsize=(14, 10))

for i, (name, sst, sss, dic, alk) in enumerate(regions):

    # ---- CLIMATOLOGY ----
    sst_m = sst.groupby('time.month').mean('time').values
    sst_s = sst.groupby('time.month').std('time').values

    sss_m = sss.groupby('time.month').mean('time').values
    sss_s = sss.groupby('time.month').std('time').values

    dic_m = dic.groupby('time.month').mean('time').values
    dic_s = dic.groupby('time.month').std('time').values

    alk_m = alk.groupby('time.month').mean('time').values
    alk_s = alk.groupby('time.month').std('time').values

    # =========================
    # LEFT: DIC & ALK
    # =========================
    ax = axes[i, 0]
    
    # DIC (left axis)
    ax.plot(x, dic_m, color=dic_clr, lw=2)
    ax.fill_between(x, dic_m-dic_s, dic_m+dic_s,
                    color=dic_clr, alpha=0.2)
    
    ax.set_ylabel(name, fontsize=11, fontweight='bold')
    ax.tick_params(axis='y', colors=dic_clr)
    ax.spines['left'].set_color(dic_clr)
    
    # ALK (right axis)
    ax2 = ax.twinx()
    ax2.plot(x, alk_m, color=alk_clr, lw=2, linestyle='--')
    ax2.fill_between(x, alk_m-alk_s, alk_m+alk_s,
                     color=alk_clr, alpha=0.2)
    
    ax2.tick_params(axis='y', colors=alk_clr)
    ax2.spines['right'].set_color(alk_clr)
    # =========================
    # RIGHT: SST & SSS
    # =========================
    ax = axes[i, 1]

    # SST (left axis)
    ax.plot(x, sst_m, color=sst_clr, lw=2)
    ax.fill_between(x, sst_m-sst_s, sst_m+sst_s,
                    color=sst_clr, alpha=0.2)
    
    ax.tick_params(axis='y', colors=sst_clr)
    ax.spines['left'].set_color(sst_clr)
    
    # SSS (right axis)
    ax2 = ax.twinx()
    ax2.plot(x, sss_m, color=sss_clr, lw=2, linestyle='--')
    ax2.fill_between(x, sss_m-sss_s, sss_m+sss_s,
                     color=sss_clr, alpha=0.2)
    
    ax2.tick_params(axis='y', colors=sss_clr)
    ax2.spines['right'].set_color(sss_clr)

# =========================
# COMMON FORMATTING
# =========================
for ax in axes.flatten():
    ax.set_xticks(x)
    ax.set_xticklabels(month_labels)
    ax.set_xlim(0, 11)
    ax.grid(alpha=0.3)

# Show x-label only at bottom row
for ax in axes[-1, :]:
    ax.set_xlabel('Month')

# =========================
# COLUMN TITLES (TOP)
# =========================
ax = axes[0,0]

ax.set_title('', fontsize=13)  # remove default title

ax.text(0.5, 1.05,
        'DIC ', color='black',
        ha='right', va='bottom',
        transform=ax.transAxes,
        fontsize=13, fontweight='bold')

ax.text(0.5, 1.05,
        '■', color=dic_clr,
        ha='left', va='bottom',
        transform=ax.transAxes,
        fontsize=13)

ax.text(0.53, 1.05,
        ' & ALK ', color='black',
        ha='left', va='bottom',
        transform=ax.transAxes,
        fontsize=13, fontweight='bold')

ax.text(0.66, 1.05,
        '■', color=alk_clr,
        ha='left', va='bottom',
        transform=ax.transAxes,
        fontsize=13)
ax = axes[0,1]

ax.set_title('', fontsize=13)

ax.text(0.5, 1.05,
        'SST ', color='black',
        ha='right', va='bottom',
        transform=ax.transAxes,
        fontsize=13, fontweight='bold')

ax.text(0.5, 1.05,
        '■', color=sst_clr,
        ha='left', va='bottom',
        transform=ax.transAxes,
        fontsize=13)

ax.text(0.53, 1.05,
        ' & SSS ', color='black',
        ha='left', va='bottom',
        transform=ax.transAxes,
        fontsize=13, fontweight='bold')

ax.text(0.66, 1.05,
        '■', color=sss_clr,
        ha='left', va='bottom',
        transform=ax.transAxes,
        fontsize=13)
# =========================
# FINAL LAYOUT
# =========================
plt.tight_layout(rect=[0.05, 0.05, 1, 0.96])
plt.show()
