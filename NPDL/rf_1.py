# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 14:38:51 2026

@author: user
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs 
import cartopy.feature as cf
import cmaps
from scipy.stats import linregress

ds = xr.open_dataset(r'/home/bobco-08/24cl05012/npdl/rf_monthly_1982_2020.nc')

ds =ds.where(ds!= 0,np.nan)

t = ds['TIME']
lat = ds['LATITUDE']
lon = ds['LONGITUDE']

rf = ds['RAINFALL']

#%%   annual climatology

rf_annual  = rf.mean(dim =('TIME')) #annual climatology
rf_av = rf_annual.values

#%% seasonal climatology function

def seasonal_climatology(da,time_dim = 'time'):
    
        djf =  da.resample(TIME='QS-DEC').mean(dim=('TIME'))
        
        djf = djf.isel(TIME=slice(1, None)) 
                
        djf = djf.sel(TIME= djf.TIME.dt.month == 12)
                
        djf = djf.mean(dim=('TIME'))

        
        mam = da.sel(TIME = da.TIME.dt.month.isin([3,4,5])).mean(dim=('TIME'))
        
        jjas = da.sel(TIME = da.TIME.dt.month.isin([6,7,8,9])).mean(dim=('TIME'))
        
        on = da.sel(TIME = da.TIME.dt.month.isin([10,11])).mean(dim=('TIME'))
        
        return { 'DJF': djf,
                  'MAM': mam,
                 'JJAS': jjas,
                  'ON': on}    

rf_season = seasonal_climatology(rf, time_dim='TIME')

rf_djf  = rf_season['DJF']
rf_mam  = rf_season['MAM']
rf_jjas = rf_season['JJAS']
rf_on   = rf_season['ON']
#%% MONTHLY grouping

rf_mon = rf.groupby('TIME.month').mean()

#%%  annual_climatology plot

fig = plt.figure(figsize=(15,10),dpi = 600)
ax = plt.axes(projection = ccrs.PlateCarree())
im = plt.contourf(lon,lat,rf_annual,cmap = cmaps.WhBlGrYeRe,levels =20,transform = ccrs.PlateCarree())
ax.coastlines(zorder = 12)
ax.add_feature(cf.LAND,color = 'grey')
gl = ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.2)

ax.set_xticks(np.arange(65,105,5), crs=ccrs.PlateCarree()) # for adding ticks 
ax.set_yticks(np.arange(5,45,5), crs=ccrs.PlateCarree())

ax.tick_params(axis='both',
               direction='out',   # makes arrow-like ticks
               length=5,
               width=1.2,
               labelsize=16)
gl.top_labels = False
gl.right_labels = False

cbar_ax = fig.add_axes([0.79, 0.15, 0.02, 0.7]) #manual colorbar
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', shrink=0.8,)
cbar.ax.tick_params(axis='both', which='major',
               direction='out',   # ticks pointing outward
               length=6,          # tick length
               width=1,labelsize=18)
cbar.set_label('Rainfall (mm)',fontsize = 18,fontweight = 'bold')
# plt.setp(ax.get_xticklabels(), fontweight='bold')
# plt.setp(ax.get_yticklabels(), fontweight='bold')
ax.set_title('Annual Climatology of Rainfall in India (1982-2020)',fontsize = 18,fontweight='bold')
plt.show()
#%% seasonal _plot

# ---- Seasons ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']

# ---- Plot configurations ----
plot_configs = [
    {
        'fields': {
            'DJF': rf_djf,
            'MAM':  rf_mam,
            'JJAS':  rf_jjas,
            'ON':   rf_on
        },
        'levels':20,
        'cmap': cmaps.WhiteBlueGreenYellowRed,
        'title': 'Rainfall in India (1982-2020)',
        'unit': 'mm'
    }
]

# ---- Main plotting loop ----
for cfg in plot_configs:

    fig, axs = plt.subplots(
        2, 2, figsize=(20,18),dpi = 600,
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': 0.03, 'hspace': 0.09})

    for i, ax in enumerate(axs.flat):

        season = seasons[i]
        data = cfg['fields'][season]

        im = ax.contourf(
            data.LONGITUDE,data.LATITUDE,data,
            levels=cfg['levels'],
            cmap=cfg['cmap'],
            transform=ccrs.PlateCarree(),
            extend='both'
        )
        
        # cs = ax.contour(
        #     data.LONGITUDE,data.LATITUDE,data,
        #     levels= 5,
        #     colors='k',
        #     linewidths=1
        # )
        
        # ax.clabel(cs, inline=True, fontsize=8, fmt='%d')

        ax.coastlines(zorder=15)
        ax.add_feature(cf.LAND, facecolor='grey')

        ax.text(
            0.85, 0.98, f"Mean ({season})",
            transform=ax.transAxes,
            fontsize= 20, zorder = 15,fontweight='bold',
            va='top', ha='right'
        )

        ax.set_xticks(np.arange(65,105,5), crs=ccrs.PlateCarree()) # for adding ticks 
        ax.set_yticks(np.arange(5,45,5), crs=ccrs.PlateCarree())

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

    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(cfg['unit'],fontsize=20,fontweight = 'bold')
    cbar.ax.tick_params(labelsize=20)

    plt.suptitle(
        f'Seasonal climatology of {cfg["title"]}',
        fontsize=21, y=0.91,fontweight = 'bold'
    )

    plt.show()
#%% monthly

fig, axs = plt.subplots(3, 4, figsize=(22,18),dpi =500,subplot_kw={'projection': ccrs.PlateCarree()})


months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']


for i, ax in enumerate(axs.flat):
    
    # Use the new levels and add extend='both'
    mon = ax.contourf(rf_mon.LONGITUDE,rf_mon.LATITUDE,rf_mon[i], 
                     levels= 30, 
                     cmap= cmaps.WhiteBlueGreenYellowRed,
                     extend='max',   
                     transform=ccrs.PlateCarree())
    
    ax.text(0.57, 0.89, months[i], 
            transform=ax.transAxes,  # Use axes coordinates (0-1)
            fontsize=20, 
            fontweight='bold',
            color='black',
            zorder = 20,
            )
    ax.coastlines(zorder=15)
    # 2. Add land and borders for better context
    ax.add_feature(cf.LAND, color='grey')

    # 3. Use ax.gridlines for cleaner label control
   # Set tick locations
    ax.set_xticks(np.arange(65,105,5), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(5,45,5), crs=ccrs.PlateCarree())
    
    # Tick style
    ax.tick_params(axis='both',
                   direction='out',
                   length=5,
                   width=1.2,
                   labelsize=20)
    
    # Show ticks only on outer plots
    ax.tick_params(labelleft=(i % 4 == 0))   # first column
    ax.tick_params(labelbottom=(i >= 8))     # bottom row
        
cbar_ax = fig.add_axes([1.01, 0.15, 0.02, 0.7])
# 2. Draw the colorbar in that new axis
cbar = fig.colorbar(mon, cax=cbar_ax,ticks = np.arange(0,300,20))
cbar.ax.tick_params(labelsize=22)
cbar.set_label(label='mm', labelpad=10, fontsize = 24,fontweight = 'bold')
# 4. Corrected main title to match your data
plt.suptitle('Monthly Mean Rainfall in India (1982–2020)', fontsize=24,y=0.97,fontweight = 'bold')
plt.tight_layout()

#%%  annual trend
annual = rf  # or rf_year if you want annual mean

time_float = ((annual.TIME - annual.TIME.isel(TIME=0)) / np.timedelta64(1, "D")).values / 365.25
time_float = time_float.astype(np.float64)
rf_float_time = annual.assign_coords(TIME=("TIME", time_float))

def linregress_1d(y, x):
    mask = np.isfinite(y) #masking all the NaN
    if mask.sum() < 2:  #regression needs more than two points minimum 
        return np.nan, np.nan
    slope, intercept, r_value, p_value, std_err = linregress(x[mask], y[mask])  # outputs
    return slope, p_value

slope, pval = xr.apply_ufunc(
    linregress_1d, rf_float_time,
    input_core_dims=[['TIME']],
    kwargs={'x': time_float},
    output_core_dims=[[], []],
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float, float]
)

rf_trend_10yr = slope * 10

signif_mask = pval < 0.05
rf_trend_signif = rf_trend_10yr.where(signif_mask)

#%%  annual trend

proj = ccrs.PlateCarree()
fig1, ax1 = plt.subplots(figsize=(14, 10),dpi = 300,subplot_kw={"projection": proj})
levels = np.linspace(-40,40,17)
map1 = ax1.contourf(lon,lat,rf_trend_10yr,cmap= cmaps.MPL_RdBu,levels = levels,transform=ccrs.PlateCarree())
sig_y, sig_x = np.where(signif_mask)  # get indices of significant points
ax1.scatter(
    lon[sig_x], lat[sig_y],
    color='k', s=2, marker='o',
    transform=ccrs.PlateCarree(),
    zorder=15,
    alpha=0.6
)
gl = ax1.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.2)

ax1.set_xticks(np.arange(65,105,5), crs=ccrs.PlateCarree()) # for adding ticks 
ax1.set_yticks(np.arange(5,45,5), crs=ccrs.PlateCarree())

ax1.tick_params(axis='both',
               direction='out',   # makes arrow-like ticks
               length=5,
               width=1.2,
               labelsize=16)

ax1.add_feature(cf.LAND, facecolor="grey")
ax1.coastlines(linewidth=0.8)
ax1.set_title("Annual rainfall trend in India in decade", fontsize=22,fontweight = 'bold')


cbar_ax = fig1.add_axes([0.83, 0.15, 0.02, 0.7])
cbar = fig.colorbar(map1, cax=cbar_ax,ticks = levels)
cbar.ax.tick_params(labelsize=18)
cbar.set_label(label="$\mathrm{mm\,decade^{-1}}$", labelpad=10, fontsize = 22,fontweight = 'bold')

plt.show()

#%% seasonal trend fuction

def seasonal(da,time_dim = 'time'):
    

        djf =  da.resample(TIME='QS-DEC').mean(dim=('TIME'))
        
        djf = djf.isel(TIME=slice(1, None)) 
                
        djf = djf.sel(TIME= djf.TIME.dt.month == 12)
        
        mam = da.sel(TIME = rf.TIME.dt.month.isin([3,4,5]))
        jjas = da.sel(TIME = rf.TIME.dt.month.isin([6,7,8,9]))
        
        on = da.sel(TIME = rf.TIME.dt.month.isin([10,11]))
        
        return { 'DJF': djf,
                  'MAM': mam,
                 'JJAS': jjas,
                  'ON': on}    

rf_season = seasonal(rf, time_dim='TIME')

rf_djf_3d  = rf_season['DJF']
rf_mam_3d  = rf_season['MAM']
rf_mam_3d = rf_mam_3d.resample(TIME='1Y').mean(dim=('TIME'))
rf_jjas_3d = rf_season['JJAS']
rf_jjas_3d = rf_jjas_3d.resample(TIME='1Y').mean(dim=('TIME'))
rf_on_3d   = rf_season['ON']
rf_on_3d = rf_on_3d.resample(TIME='1Y').mean(dim=('TIME'))
#%% seasonal trends

djf_float = ((rf_djf_3d.TIME - rf_djf_3d.TIME.isel(TIME=0)) / np.timedelta64(1, "D")).values / 365.25
djf_float = djf_float.astype(np.float64)
rf_djf_time = rf_djf_3d.assign_coords(TIME=("TIME", djf_float))
# ===================== DJF =====================

djf_slope, djf_pval = xr.apply_ufunc(
    linregress_1d, rf_djf_time,
    input_core_dims=[['TIME']],
    kwargs={'x': djf_float},
    output_core_dims=[[], []],
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float, float]
)

djf_trend_10yr = djf_slope * 10

djf_signif_mask = djf_pval < 0.05
djf_trend_signif = djf_trend_10yr.where(signif_mask)
#%%
mam_float = ((rf_mam_3d.TIME - rf_mam_3d.TIME.isel(TIME=0)) / np.timedelta64(1, "D")).values / 365.25
mam_float = mam_float.astype(np.float64)
rf_mam_time = rf_mam_3d.assign_coords(TIME=("TIME", mam_float))

mam_slope, mam_pval = xr.apply_ufunc(
    linregress_1d, rf_mam_time,
    input_core_dims=[['TIME']],
    kwargs={'x': mam_float},
    output_core_dims=[[], []],
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float, float]
)

mam_trend_10yr = mam_slope * 10

mam_signif_mask = mam_pval < 0.05
mam_trend_signif = mam_trend_10yr.where(signif_mask)
#%%

jjas_float = ((rf_jjas_3d.TIME - rf_jjas_3d.TIME.isel(TIME=0)) / np.timedelta64(1, "D")).values / 365.25
jjas_float = jjas_float.astype(np.float64)
rf_jjas_time = rf_jjas_3d.assign_coords(TIME=("TIME", jjas_float))

jjas_slope, jjas_pval = xr.apply_ufunc(
    linregress_1d, rf_jjas_time,
    input_core_dims=[['TIME']],
    kwargs={'x': jjas_float},
    output_core_dims=[[], []],
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float, float]
)

jjas_trend_10yr = jjas_slope * 10

jjas_signif_mask = jjas_pval < 0.05
jjas_trend_signif = jjas_trend_10yr.where(signif_mask)
#%%
on_float = ((rf_on_3d.TIME - rf_on_3d.TIME.isel(TIME=0)) / np.timedelta64(1, "D")).values / 365.25
on_float = on_float.astype(np.float64)
rf_on_time = rf_on_3d.assign_coords(TIME=("TIME", on_float))

on_slope, on_pval = xr.apply_ufunc(
    linregress_1d, rf_on_time,
    input_core_dims=[['TIME']],
    kwargs={'x': on_float},
    output_core_dims=[[], []],
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float, float]
)

on_trend_10yr = on_slope * 10

on_signif_mask = on_pval < 0.05
on_trend_signif = on_trend_10yr.where(signif_mask)




#%% seasonal trend plots

# ---- Seasons ----
seasons = ['DJF','MAM','JJAS','ON']

# ---- Data and significance masks ----
season_data = [djf_trend_10yr, mam_trend_10yr, jjas_trend_10yr, on_trend_10yr]
season_sig  = [djf_signif_mask, mam_signif_mask, jjas_signif_mask, on_signif_mask]

levels = np.linspace(-40,40,11)

fig, axs = plt.subplots(
    2,2, figsize=(20,18),dpi = 500,
    subplot_kw={'projection': ccrs.PlateCarree()},
    gridspec_kw={'wspace':-0.10,'hspace':0.04}
)

for i, ax in enumerate(axs.flat):

    data = season_data[i]
    signif_mask = season_sig[i]
    season = seasons[i]

    im = ax.contourf(
        data.LONGITUDE, data.LATITUDE, data,
        levels=levels,
        cmap=cmaps.MPL_RdBu,
        transform=ccrs.PlateCarree(),
        extend='both'
    )

    # ---- Significance dots ----
    sig_y, sig_x = np.where(signif_mask)

    ax.scatter(
        data.LONGITUDE[sig_x],
        data.LATITUDE[sig_y],
        color='k',
        s=2,
        marker='o',
        transform=ccrs.PlateCarree(),
        zorder=15,
        alpha=0.6
    )

    ax.coastlines(zorder=15)
    ax.add_feature(cf.LAND, facecolor='grey')

    ax.text(
        0.85,0.98,season,
        transform=ax.transAxes,
        fontsize=12,
        fontweight='bold',
        va='top',ha='right'
    )
    ax.set_xticks(np.arange(65,105,5), crs=ccrs.PlateCarree()) # for adding ticks 
    ax.set_yticks(np.arange(5,45,5), crs=ccrs.PlateCarree())
    
    ax.tick_params(axis='both',
                   direction='out',   # makes arrow-like ticks
                   length=5,
                   width=1.2,
                   labelsize=18)
    
      # Hide labels depending on subplot index
    if i in [0, 1]:          # top row
        ax.tick_params(labelbottom=False)
    
    if i in [1, 3]:          # right column
        ax.tick_params(labelleft=False)
   

# ---- Colorbar ----
cbar_ax = fig.add_axes([0.90,0.15,0.02,0.7])
cbar = fig.colorbar(im, cax=cbar_ax, ticks=levels)

cbar.set_label('mm decade$^{-1}$', fontsize=22,fontweight = 'bold')
cbar.ax.tick_params(labelsize=18)

plt.suptitle('Seasonal rainfall trend in India in decade', fontsize=22, y=0.92,fontweight ='bold')

plt.show()
    
#%% correlation rf vs iod

ds_iod = xr.open_dataset('/home/bobco-08/24cl05012/npdl/dmi.had.long.nc')    

iod = ds_iod['value']
iod =iod.sel(time =slice('1982-01-01','2020-12-01'))
iod = iod.rename({'time': 'TIME'})

corr_rf_iod = xr.corr(rf,iod,dim ='TIME')

#%%

proj = ccrs.PlateCarree()
fig1, ax1 = plt.subplots(figsize=(14,12),dpi=300,subplot_kw={"projection": proj})
levels = np.linspace(-0.20,0.20,11)
map1 = ax1.contourf(lon,lat,corr_rf_iod,cmap= cmaps.MPL_RdBu_r,levels = levels,transform=ccrs.PlateCarree())
ax1.scatter(
    lon[sig_x], lat[sig_y],
    color='k', s=2, marker='o',
    transform=ccrs.PlateCarree(),
    zorder=15,
    alpha=0.8
)
ax1.set_xticks(np.arange(65,105,5), crs=ccrs.PlateCarree()) # for adding ticks 
ax1.set_yticks(np.arange(5,45,5), crs=ccrs.PlateCarree())

ax1.tick_params(axis='both',
               direction='out',   # makes arrow-like ticks
               length=5,
               width=1.2,
               labelsize=18)

ax1.add_feature(cf.LAND, facecolor="grey")
ax1.coastlines(linewidth=0.8)
ax1.set_title("Correlation of Rainfall vs IOD in Indian Region (1982-2020)", fontsize=20,fontweight ='bold')
cb1 = plt.colorbar(map1, ax=ax1, orientation="vertical", fraction=0.045,pad=0.05,aspect=40, ticks =levels)
cb1.set_label("Correlation", fontsize=18,fontweight ='bold')
cb1.ax.tick_params(labelsize=18)
plt.show()

#%%

ds_enso = xr.open_dataset(r'/home/bobco-08/24cl05012/npdl/nino34.long.anom.nc')    

enso = ds_enso['value']
enso = enso.sel(time =slice('1982-01-01','2020-12-01'))
enso = enso.rename({'time': 'TIME'})

corr_rf_enso = xr.corr(rf,enso,dim ='TIME')
#%%

proj = ccrs.PlateCarree()
fig1, ax1 = plt.subplots(figsize=(14,12),dpi = 300,subplot_kw={"projection": proj})
levels = np.linspace(-0.20,0.20,11)
map1 = ax1.contourf(lon,lat,corr_rf_enso,cmap= cmaps.MPL_RdBu_r,levels = levels,transform=ccrs.PlateCarree())
ax1.scatter(
    lon[sig_x], lat[sig_y],
    color='k', s=2, marker='o',
    transform=ccrs.PlateCarree(),
    zorder=15,
    alpha=0.6
)
ax1.set_xticks(np.arange(65,105,5), crs=ccrs.PlateCarree()) # for adding ticks 
ax1.set_yticks(np.arange(5,45,5), crs=ccrs.PlateCarree())

ax1.tick_params(axis='both',
               direction='out',   # makes arrow-like ticks
               length=5,
               width=1.2,
               labelsize=18)

ax1.add_feature(cf.LAND, facecolor="grey")
ax1.coastlines(linewidth=0.8)
ax1.set_title("Correlation of Rainfall vs ENSO in Indian Region (1982-2020)", fontsize=20,fontweight ='bold')
cb1 = plt.colorbar(map1, ax=ax1, orientation="vertical", fraction=0.045,pad=0.05,aspect=40, ticks =levels)
cb1.set_label("Correlation", fontsize=18,fontweight ='bold')
cb1.ax.tick_params(labelsize=18)
plt.show()
#%%
import ruptures as rpt
import matplotlib.pyplot as plt

rf_t = rf.mean(dim=('LATITUDE','LONGITUDE'))

signal_ls = rf_t.values.astype(float).reshape(-1,1)
time_ls = rf_t['TIME'].values

# Binary segmentation with l2 model
algo = rpt.Binseg(model="l2", min_size=48).fit(signal_ls)

# Force exactly 3 change points
result_gs = algo.predict(n_bkps=3)
gs_yrs = np.array(result_gs)



for i in range(len(gs_yrs)-1):
    regime_years  = time_ls[gs_yrs[i]]
    print(regime_years)
    
plt.figure(figsize=(12,5))
plt.plot(time_ls, signal_ls[:,0], 'b', lw=1.8, label="rainfall")

for bp in result_gs[:-1]:
    plt.axvline(time_ls[bp], color='r', ls='--', lw=1.8, alpha=1)
plt.title(f"Change Point Detection in annual rainfall in India (5°N – 38°N,65°E – 100°E)", fontsize=16)

plt.xlabel("Year")
plt.ylabel("Rainfall (mm)")

plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
#%%
regime_1 = rf.sel(TIME = slice('1982-01-31','1987-06-01')).mean(dim=('TIME'))
regime_2 = rf.sel(TIME = slice('1987-06-01','1995-10-01')).mean(dim=('TIME'))
regime_3 = rf.sel(TIME = slice('1995-10-01','2005-05-01')).mean(dim=('TIME'))
regime_4 = rf.sel(TIME = slice('2005-05-01','2020-12-01')).mean(dim=('TIME'))

#%%

regimes = [regime_1, regime_2, regime_3, regime_4]
titles = ['1982–1987', '1987–1995', '1995–2005', '2005–2020']

fig, ax = plt.subplots(
    2, 2, figsize=(14,10),dpi =500,
    subplot_kw={'projection': ccrs.PlateCarree()} ,gridspec_kw={'wspace': -0.28, 'hspace': 0.20}
 
)

for i, r in enumerate(regimes):
    a = ax.flat[i]
    
    im = a.contourf(
        rf.LONGITUDE, rf.LATITUDE, r,
        levels=35,
        cmap=cmaps.WhiteBlueGreenYellowRed,
        transform=ccrs.PlateCarree(),
        extend='both'
    )
    
    # Coastlines and land
    a.coastlines()
    a.add_feature(cf.LAND, facecolor='lightgrey')

    a.set_xticks(np.arange(65,105,5))
    a.set_yticks(np.arange(5,45,5))

    a.tick_params(axis='both', direction='out', length=5, width=1.2, labelsize=12)
    
    # Subplot title
    a.set_title(titles[i], fontsize=14)

# ---- Manual colorbar ----
cbar_ax = fig.add_axes([0.82,0.15,0.02,0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label('mm', fontsize=18)
cbar.ax.tick_params(labelsize=18)

# ---- Overall title ----
plt.suptitle('Annual Climatology of Rainfall in different regimes', fontsize=18, fontweight='bold', y=0.94)

plt.tight_layout()
plt.show()

#%%

regime_1 = rf.sel(TIME = slice('1982-01-31','1987-06-01'))
regime_2 = rf.sel(TIME = slice('1987-06-01','1995-10-01'))
regime_3 = rf.sel(TIME = slice('1995-10-01','2005-05-01'))
regime_4 = rf.sel(TIME = slice('2005-05-01','2020-12-01'))

# ---- List of regimes ----
regimes = [regime_1, regime_2, regime_3, regime_4]
regime_titles = ['1982–1987', '1987–1995', '1995–2005', '2005–2020']

# ---- Seasons ----
seasons = ['DJF', 'MAM', 'JJAS', 'ON']
 
# ---- Loop over each regime ----
for idx, regime in enumerate(regimes):
    
    # Compute seasonal climatology
    rf_season = seasonal_climatology(regime, time_dim='TIME')
    
    # Create 2x2 subplot for this regime
    fig, ax = plt.subplots(
        2, 2, figsize=(10,10),dpi=500,
        subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Plot each season
    for i, season in enumerate(seasons):
        a = ax.flat[i]
        data = rf_season[season]
        
        # Contourf
        im = a.contourf(
            data.LONGITUDE, data.LATITUDE, data,
            levels=levels,
            cmap= cmaps.WhiteBlueGreenYellowRed,
            transform=ccrs.PlateCarree(),
            extend='both'
        )
        
      
        a.coastlines()
        a.add_feature(cf.LAND, facecolor='lightgrey')
        
        # Custom ticks
        a.set_xticks(np.arange(65,105,5))
        a.set_yticks(np.arange(5,45,5))
        
        if i in [0, 1]:          # top row
            a.tick_params(labelbottom=False)
        
        if i in [1, 3]:          # right column
            a.tick_params(labelleft=False)
        # Arrow-like ticks
        a.tick_params(axis='both', direction='out', length=4, width=1.2, labelsize=10)
        
        # Subplot title
        a.set_title(seasons[i], fontsize=14, fontweight='bold')
    
    # ---- Manual colorbar ----
    cbar_ax = fig.add_axes([1,0.15,0.02,0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Rainfall (mm)', fontsize=14, fontweight='bold')
    cbar.ax.tick_params(labelsize=12)
    
    # ---- Overall figure title ----
    plt.suptitle(f'Seasonal Rainfall Climatology in Regime {regime_titles[idx]}',
                 fontsize=20, fontweight='bold', y=0.99)
    
    plt.tight_layout()
    plt.show()