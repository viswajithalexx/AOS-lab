#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 15:54:35 2026

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


file = '/home/bobco-08/24cl05012/CO2/data/data_1/cmems/bio_var/cmems_obs-mob_glo_bgc-car_my_irr-i_1775027250539_30E-119.88E_29.88S-29.88N.nc'

data = xr.open_dataset(file)
data = data.rename({'latitude':'lat', 'longitude': 'lon'})
lat = data['lat']
lon =data['lon']
time = data['time']

data2 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/cmems/phy_var/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_multi-vars_35.00E-120.00E_30.00S-30.00N_0.51m_1993-01-01-2024-12-01.nc')
data2 = data2.rename({'latitude':'lat', 'longitude': 'lon'})
sst = data2['thetao_oras'].isel(depth=0)
sss = data2['so_oras'].isel(depth=0)


#%% nwio

dic_nwio = data['tco2'].sel(lat = slice(5,22.5),lon = slice(45,65),time = slice('1994-01-01','2024-12-01') )
alk_nwio = data['talk'].sel(lat = slice(5,22.5),lon = slice(45,65),time = slice('1994-01-01','2024-12-01') )
pco2_nwio = data['spco2'].sel(lat = slice(5,22.5),lon = slice(45,65),time = slice('1994-01-01','2024-12-01') )
sst_nwio = sst.sel(lat = slice(5,22.5),lon = slice(45,65),time = slice('1994-01-01','2024-12-01') )
sss_nwio = sss.sel(lat = slice(5,22.5),lon = slice(45,65),time = slice('1994-01-01','2024-12-01'))

#%% eio

dic_eio = data['tco2'].sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01') )
alk_eio = data['talk'].sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01') )
pco2_eio = data['spco2'].sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01') )
sst_eio = sst.sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01') )
sss_eio = sss.sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01'))

#%% nas

dic_nas = data['tco2'].sel(lat = slice(22.5,28),lon = slice(56,70),time = slice('1994-01-01','2024-12-01') )
alk_nas = data['talk'].sel(lat = slice(22.5,28),lon = slice(56,70),time = slice('1994-01-01','2024-12-01') )
pco2_nas = data['spco2'].sel(lat = slice(22.5,28),lon = slice(56,70),time = slice('1994-01-01','2024-12-01') )
sst_nas = sst.sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01') )
sss_nas = sss.sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01') )

#%% esio

dic_esio = data['tco2'].sel(lat = slice(-6.6,8),lon = slice(92,109),time = slice('1994-01-01','2024-12-01') )
alk_esio = data['talk'].sel(lat = slice(-6.6,8),lon = slice(92,109),time = slice('1994-01-01','2024-12-01') )
pco2_esio = data['spco2'].sel(lat = slice(-6.6,8),lon = slice(92,109),time = slice('1994-01-01','2024-12-01') )
sst_esio = sst.sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01') )
sss_esio = sss.sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01') )

#%% scatter plot -- nwio


# Spatial mean time series
pco2_nwio_t = pco2_nwio.mean(dim=('lat','lon'))

dic_nwio_t = dic_nwio.mean(dim=('lat','lon'))
alk_nwio_t = alk_nwio.mean(dim=('lat','lon'))
sst_nwio_t = sst_nwio.mean(dim=('lat','lon'))
sss_nwio_t = sss_nwio.mean(dim=('lat','lon'))

# Plot
# Correlations using xarray
r_dic_nwio = xr.corr(dic_nwio_t, pco2_nwio_t, dim='time').values
r_alk_nwio = xr.corr(alk_nwio_t, pco2_nwio_t, dim='time').values
r_sst_nwio = xr.corr(sst_nwio_t, pco2_nwio_t, dim='time').values
r_sss_nwio = xr.corr(sss_nwio_t, pco2_nwio_t, dim='time').values

fig, ax = plt.subplots(2,2, figsize=(18,14))

fig.suptitle('Northwestern Indian Ocean (NWIO)', fontsize=16, fontweight='bold', y=0.98)

# DIC
ax[0,0].scatter(dic_nwio_t, pco2_nwio_t, s=28)
ax[0,0].set_xlabel('DIC (µmol kg$^{-1}$)', fontsize=20)
ax[0,0].set_ylabel('pCO$_2$ (µatm)', fontsize=20)
ax[0,0].set_title(f'DIC vs pCO$_2$ (r = {r_dic_nwio:.2f})', fontsize=20)
ax[0,0].tick_params(axis='both', labelsize=16)

# ALK
ax[0,1].scatter(alk_nwio_t, pco2_nwio_t, s=28)
ax[0,1].set_xlabel('ALK (µmol kg$^{-1}$)', fontsize=20)
ax[0,1].set_ylabel('pCO$_2$ (µatm)', fontsize=20)
ax[0,1].set_title(f'ALK vs pCO$_2$ (r = {r_alk_nwio:.2f})', fontsize=20)
ax[0,1].tick_params(axis='both', labelsize=16)

# SST
ax[1,0].scatter(sst_nwio_t, pco2_nwio_t, s=28)
ax[1,0].set_xlabel('SST (°C)', fontsize=20)
ax[1,0].set_ylabel('pCO$_2$ (µatm)', fontsize=20)
ax[1,0].set_title(f'SST vs pCO$_2$ (r = {r_sst_nwio:.2f})', fontsize=20)
ax[1,0].tick_params(axis='both', labelsize=16)

# SSS
ax[1,1].scatter(sss_nwio_t, pco2_nwio_t, s=28)
ax[1,1].set_xlabel('SSS (psu)', fontsize=20)
ax[1,1].set_ylabel('pCO$_2$ (µatm)', fontsize=20)
ax[1,1].set_title(f'SSS vs pCO$_2$ (r = {r_sss_nwio:.2f})', fontsize=20)
ax[1,1].tick_params(axis='both', labelsize=16)

plt.tight_layout()
plt.show()

#%% scatter plot -- eio

# Spatial mean time series
pco2_eio_t = pco2_eio.mean(dim=('lat','lon'))
dic_eio_t = dic_eio.mean(dim=('lat','lon'))
alk_eio_t = alk_eio.mean(dim=('lat','lon'))
sst_eio_t = sst_eio.mean(dim=('lat','lon'))
sss_eio_t = sss_eio.mean(dim=('lat','lon'))

# Correlations
r_dic_eio = xr.corr(dic_eio_t, pco2_eio_t, dim='time').values
r_alk_eio = xr.corr(alk_eio_t, pco2_eio_t, dim='time').values
r_sst_eio = xr.corr(sst_eio_t, pco2_eio_t, dim='time').values
r_sss_eio = xr.corr(sss_eio_t, pco2_eio_t, dim='time').values

# Plot
fig, ax = plt.subplots(2,2, figsize=(10,8))

# Common title
fig.suptitle('Equatorial Indian Ocean (EIO)', fontsize=16, fontweight='bold',y=0.98)

# DIC
ax[0,0].scatter(dic_eio_t, pco2_eio_t, s=10)
ax[0,0].set_xlabel('DIC (µmol kg$^{-1}$)')
ax[0,0].set_ylabel('pCO$_2$ (µatm)')
ax[0,0].set_title(f'DIC vs pCO$_2$ (r = {r_dic_eio:.2f})')

# ALK
ax[0,1].scatter(alk_eio_t, pco2_eio_t, s=10)
ax[0,1].set_xlabel('ALK (µmol kg$^{-1}$)')
ax[0,1].set_ylabel('pCO$_2$ (µatm)')
ax[0,1].set_title(f'ALK vs pCO$_2$ (r = {r_alk_eio:.2f})')

# SST
ax[1,0].scatter(sst_eio_t, pco2_eio_t, s=10)
ax[1,0].set_xlabel('SST (°C)')
ax[1,0].set_ylabel('pCO$_2$ (µatm)')
ax[1,0].set_title(f'SST vs pCO$_2$ (r = {r_sst_eio:.2f})')

# SSS
ax[1,1].scatter(sss_eio_t, pco2_eio_t, s=10)
ax[1,1].set_xlabel('SSS (psu)')
ax[1,1].set_ylabel('pCO$_2$ (µatm)')
ax[1,1].set_title(f'SSS vs pCO$_2$ (r = {r_sss_eio:.2f})')

plt.tight_layout()
plt.show()

#%% scatter plot - nas

# Spatial mean time series
pco2_nas_t = pco2_nas.mean(dim=('lat','lon'))

dic_nas_t = dic_nas.mean(dim=('lat','lon'))
alk_nas_t = alk_nas.mean(dim=('lat','lon'))
sst_nas_t = sst_nas.mean(dim=('lat','lon'))
sss_nas_t = sss_nas.mean(dim=('lat','lon'))

# Correlations
r_dic_nas = xr.corr(dic_nas_t, pco2_nas_t, dim='time').values
r_alk_nas = xr.corr(alk_nas_t, pco2_nas_t, dim='time').values
r_sst_nas = xr.corr(sst_nas_t, pco2_nas_t, dim='time').values
r_sss_nas = xr.corr(sss_nas_t, pco2_nas_t, dim='time').values

# Plot
fig, ax = plt.subplots(2,2, figsize=(10,8))

# Common title
fig.suptitle('Northern Arabian Sea (NAS)', fontsize=16, fontweight='bold',y=0.95)

# DIC
ax[0,0].scatter(dic_nas_t, pco2_nas_t, s=10)
ax[0,0].set_xlabel('DIC (µmol kg$^{-1}$)')
ax[0,0].set_ylabel('pCO$_2$ (µatm)')
ax[0,0].set_title(f'DIC vs pCO$_2$ (r = {r_dic_nas:.2f})')

# ALK
ax[0,1].scatter(alk_nas_t, pco2_nas_t, s=10)
ax[0,1].set_xlabel('ALK (µmol kg$^{-1}$)')
ax[0,1].set_ylabel('pCO$_2$ (µatm)')
ax[0,1].set_title(f'ALK vs pCO$_2$ (r = {r_alk_nas:.2f})')

# SST
ax[1,0].scatter(sst_nas_t, pco2_nas_t, s=10)
ax[1,0].set_xlabel('SST (°C)')
ax[1,0].set_ylabel('pCO$_2$ (µatm)')
ax[1,0].set_title(f'SST vs pCO$_2$ (r = {r_sst_nas:.2f})')

# SSS
ax[1,1].scatter(sss_nas_t, pco2_nas_t, s=10)
ax[1,1].set_xlabel('SSS (psu)')
ax[1,1].set_ylabel('pCO$_2$ (µatm)')
ax[1,1].set_title(f'SSS vs pCO$_2$ (r = {r_sss_nas:.2f})')

plt.tight_layout(rect=[0,0,1,0.95])
plt.show()

#%% scatter plots - esio

# Spatial mean time series
pco2_esio_t = pco2_esio.mean(dim=('lat','lon'))

dic_esio_t = dic_esio.mean(dim=('lat','lon'))
alk_esio_t = alk_esio.mean(dim=('lat','lon'))
sst_esio_t = sst_esio.mean(dim=('lat','lon'))
sss_esio_t = sss_esio.mean(dim=('lat','lon'))


# Correlations
r_dic_esio = xr.corr(dic_esio_t, pco2_esio_t, dim='time').values
r_alk_esio = xr.corr(alk_esio_t, pco2_esio_t, dim='time').values
r_sst_esio = xr.corr(sst_esio_t, pco2_esio_t, dim='time').values
r_sss_esio = xr.corr(sss_esio_t, pco2_esio_t, dim='time').values


# Plot
fig, ax = plt.subplots(2,2, figsize=(10,8))

fig.suptitle('Eastern Indian Ocean (ESIO)', fontsize=16, fontweight='bold',y = 0.95)

# DIC
ax[0,0].scatter(dic_esio_t, pco2_esio_t, s=10)
ax[0,0].set_xlabel('DIC (µmol kg$^{-1}$)')
ax[0,0].set_ylabel('pCO$_2$ (µatm)')
ax[0,0].set_title(f'DIC vs pCO$_2$ (r = {r_dic_esio:.2f})')

# ALK
ax[0,1].scatter(alk_esio_t, pco2_esio_t, s=10)
ax[0,1].set_xlabel('ALK (µmol kg$^{-1}$)')
ax[0,1].set_ylabel('pCO$_2$ (µatm)')
ax[0,1].set_title(f'ALK vs pCO$_2$ (r = {r_alk_esio:.2f})')

# SST
ax[1,0].scatter(sst_esio_t, pco2_esio_t, s=10)
ax[1,0].set_xlabel('SST (°C)')
ax[1,0].set_ylabel('pCO$_2$ (µatm)')
ax[1,0].set_title(f'SST vs pCO$_2$ (r = {r_sst_esio:.2f})')

# SSS
ax[1,1].scatter(sss_esio_t, pco2_esio_t, s=10)
ax[1,1].set_xlabel('SSS (psu)')
ax[1,1].set_ylabel('pCO$_2$ (µatm)')
ax[1,1].set_title(f'SSS vs pCO$_2$ (r = {r_sss_esio:.2f})')
plt.tight_layout(rect=[0,0,1,0.95])
plt.show()

#%% scatter + density - esio


from scipy.stats import gaussian_kde


fig, ax = plt.subplots(2,2, figsize=(15,13))

fig.suptitle('Eastern South Indian Ocean (ESIO)', fontsize=20, fontweight='bold', y=0.95)

def density_scatter(x, y, axis):

    x = x.values
    y = y.values
    
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    # Normalize density to 0–100 like the paper figure
    z = 100 * (z - z.min()) / (z.max() - z.min())

    sc = axis.scatter(x, y, c=z, cmap='viridis', s=50, edgecolor='none')

    return sc

# DIC
sc = density_scatter(dic_esio_t, pco2_esio_t, ax[0,0])
ax[0,0].set_xlabel('DIC (µmol kg$^{-1}$)',fontsize = 20)
ax[0,0].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[0,0].set_title(f'(a) r = {r_dic_esio:.2f}',fontsize = 20)

# ALK
sc = density_scatter(alk_esio_t, pco2_esio_t, ax[0,1])
ax[0,1].set_xlabel('ALK (µmol kg$^{-1}$)',fontsize = 20)
ax[0,1].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[0,1].set_title(f'(b) r = {r_alk_esio:.2f}',fontsize = 20)

# SST
sc = density_scatter(sst_esio_t, pco2_esio_t, ax[1,0])
ax[1,0].set_xlabel('SST (°C)',fontsize = 20)
ax[1,0].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[1,0].set_title(f'(c) r = {r_sst_nas:.2f}',fontsize = 20)

# SSS
sc = density_scatter(sss_esio_t, pco2_esio_t, ax[1,1])
ax[1,1].set_xlabel('SSS (psu)',fontsize = 20)
ax[1,1].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[1,1].set_title(f'(d) r = {r_sss_esio:.2f}',fontsize = 20)

# ---- Adjust x and y tick size
for a in ax.ravel():
    a.tick_params(axis='both', labelsize=18)
    a.tick_params(direction='in', length=5, width=1)



# Common colorbar
cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])

cbar = fig.colorbar(sc, cax=cbar_ax)
cbar.set_label('Density',fontsize = 18)
cbar.ax.tick_params(labelsize=18)

plt.tight_layout(rect=[0,0.05,1,0.95])
plt.show()

#%% scatter + denasity - nwio
fig, ax = plt.subplots(2,2, figsize=(15,13))

fig.suptitle('Northwestern Indian Ocean (NWIO)', fontsize=20, fontweight='bold', y=0.95)
# DIC
sc = density_scatter(dic_nwio_t, pco2_nwio_t, ax[0,0])
ax[0,0].set_xlabel('DIC (µmol kg$^{-1}$)',fontsize = 20)
ax[0,0].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[0,0].set_title(f'(a) r = {r_dic_nwio:.2f}',fontsize = 20)

# ALK
sc = density_scatter(alk_nwio_t, pco2_nwio_t, ax[0,1])
ax[0,1].set_xlabel('ALK (µmol kg$^{-1}$)',fontsize = 20)
ax[0,1].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[0,1].set_title(f'(b) r = {r_alk_nwio:.2f}',fontsize = 20)

# SST
sc = density_scatter(sst_nwio_t, pco2_nwio_t, ax[1,0])
ax[1,0].set_xlabel('SST (°C)',fontsize = 20)
ax[1,0].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[1,0].set_title(f'(c) r = {r_sst_nwio:.2f}',fontsize = 20)

# SSS
sc = density_scatter(sss_nwio_t, pco2_nwio_t, ax[1,1])
ax[1,1].set_xlabel('SSS (psu)',fontsize = 20)
ax[1,1].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[1,1].set_title(f'(d) r = {r_sss_nwio:.2f}',fontsize = 20)

for a in ax.ravel():
    a.tick_params(axis='both', labelsize=18)
    a.tick_params(direction='in', length=5, width=1)

# Common colorbar
cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])

cbar = fig.colorbar(sc, cax=cbar_ax)
cbar.set_label('Density',fontsize = 18)
cbar.ax.tick_params(labelsize=18)

plt.tight_layout(rect=[0,0.05,1,0.95])
plt.show()

#%% scatter + density - eio
# DIC

fig, ax = plt.subplots(2,2, figsize=(15,13))

fig.suptitle('Equatorial Indian Ocean (EIO)', fontsize=20, fontweight='bold', y=0.95)

sc = density_scatter(dic_eio_t, pco2_eio_t, ax[0,0])
ax[0,0].set_xlabel('DIC (µmol kg$^{-1}$)',fontsize = 20)
ax[0,0].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[0,0].set_title(f'(a) r = {r_dic_eio:.2f}',fontsize = 20)

# ALK
sc = density_scatter(alk_eio_t, pco2_eio_t, ax[0,1])
ax[0,1].set_xlabel('ALK (µmol kg$^{-1}$)',fontsize = 20)
ax[0,1].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[0,1].set_title(f'(b) r = {r_alk_eio:.2f}',fontsize = 20)

# SST
sc = density_scatter(sst_eio_t, pco2_eio_t, ax[1,0])
ax[1,0].set_xlabel('SST (°C)',fontsize = 20)
ax[1,0].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[1,0].set_title(f'(c) r = {r_sst_eio:.2f}',fontsize = 20)

# SSS
sc = density_scatter(sss_eio_t, pco2_eio_t, ax[1,1])
ax[1,1].set_xlabel('SSS (psu)',fontsize = 20)
ax[1,1].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[1,1].set_title(f'(d) r = {r_sss_eio:.2f}',fontsize = 20)

for a in ax.ravel():
    a.tick_params(axis='both', labelsize=18)
    a.tick_params(direction='in', length=5, width=1)

# Common colorbar
cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])

cbar = fig.colorbar(sc, cax=cbar_ax)
cbar.set_label('Density',fontsize = 18)
cbar.ax.tick_params(labelsize=18)

plt.tight_layout(rect=[0,0.05,1,0.95])
plt.show()

#%% scatter + density - nas


fig, ax = plt.subplots(2,2, figsize=(15,13))

fig.suptitle('Northern Arabian Sea (NAS)', fontsize=20, fontweight='bold', y=0.95)


# DIC
sc = density_scatter(dic_nas_t, pco2_nas_t, ax[0,0])
ax[0,0].set_xlabel('DIC (µmol kg$^{-1}$)',fontsize = 20)
ax[0,0].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[0,0].set_title(f'(a) r = {r_dic_nas:.2f}',fontsize = 20)

# ALK
sc = density_scatter(alk_nas_t, pco2_nas_t, ax[0,1])
ax[0,1].set_xlabel('ALK (µmol kg$^{-1}$)',fontsize = 20)
ax[0,1].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[0,1].set_title(f'(b) r = {r_alk_nas:.2f}',fontsize = 20)

# SST
sc = density_scatter(sst_nas_t, pco2_nas_t, ax[1,0])
ax[1,0].set_xlabel('SST (°C)',fontsize = 20)
ax[1,0].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[1,0].set_title(f'(c) r = {r_sst_nas:.2f}',fontsize = 20)

# SSS
sc = density_scatter(sss_nas_t, pco2_nas_t, ax[1,1])
ax[1,1].set_xlabel('SSS (psu)',fontsize = 20)
ax[1,1].set_ylabel('pCO$_2$ (µatm)',fontsize = 20)
ax[1,1].set_title(f'(d) r = {r_sss_nas:.2f}',fontsize = 20)

for a in ax.ravel():
    a.tick_params(axis='both', labelsize=18)
    a.tick_params(direction='in', length=5, width=1)

# Common colorbar
cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])

cbar = fig.colorbar(sc, cax=cbar_ax)
cbar.set_label('Density',fontsize = 18)
cbar.ax.tick_params(labelsize=18)

plt.tight_layout(rect=[0,0.05,1,0.95])
plt.show()
#%%
import matplotlib.pyplot as plt

# -------- Organize data --------
regions = [
    ("NWIO", dic_nwio_t, alk_nwio_t, sst_nwio_t, sss_nwio_t, pco2_nwio_t,
     r_dic_nwio, r_alk_nwio, r_sst_nwio, r_sss_nwio),

    ("EIO", dic_eio_t, alk_eio_t, sst_eio_t, sss_eio_t, pco2_eio_t,
     r_dic_eio, r_alk_eio, r_sst_eio, r_sss_eio),

    ("NAS", dic_nas_t, alk_nas_t, sst_nas_t, sss_nas_t, pco2_nas_t,
     r_dic_nas, r_alk_nas, r_sst_nas, r_sss_nas),

    ("ESIO", dic_esio_t, alk_esio_t, sst_esio_t, sss_esio_t, pco2_esio_t,
     r_dic_esio, r_alk_esio, r_sst_esio, r_sss_esio)
]

xlabels = [
    'DIC (µmol kg$^{-1}$)',
    'ALK (µmol kg$^{-1}$)',
    'SST (°C)',
    'SSS (psu)'
]

# -------- Create figure --------
fig, ax = plt.subplots(4, 4, figsize=(32,30))

# -------- Loop --------
for i, (name, dic, alk, sst, sss, pco2, r_dic, r_alk, r_sst, r_sss) in enumerate(regions):

    variables = [dic, alk, sst, sss]
    correlations = [r_dic, r_alk, r_sst, r_sss]

    for j in range(4):

        ax[i, j].scatter(variables[j], pco2, s=20)

        # X labels only for bottom row
        if i == 3:
            ax[i, j].set_xlabel(xlabels[j], fontsize=14)

        # Y label only for first column
        if j == 0:
            ax[i, j].set_ylabel(f'{name}\npCO$_2$ (µatm)', fontsize=24)

        # Correlation title
        ax[i, j].set_title(f'r = {correlations[j]:.2f}', fontsize=24)

        # Tick style
        ax[i, j].tick_params(axis='both', labelsize=24, direction='in')

# -------- Layout --------
plt.tight_layout(rect=[0,0,1,0.92])
plt.show()