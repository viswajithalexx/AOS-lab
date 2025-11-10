#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 13:17:53 2025

@author: bobco-08
"""

import xarray as xr

file0 = '/home/bobco-08/24cl05012/CO2/data/pco2_climatology_taka.nc'
file1 = "/home/bobco-08/24cl05012/CO2/data/oc_v2024E.pCO2.nc"
file2 = '/home/bobco-08/24cl05012/CO2/data/pCO2-Corrected_INCOIS-BIO-ROMS_v2.nc'

data0 = xr.open_dataset(file0)
data1 = xr.open_dataset(file1)
data2 = xr.open_dataset(file2)

lat0 = data0['lat']
lon0 = data0['lon']
mon = data0['month']
pco2_0 = data0['pco2_sw']

lat1 = data1['lat']
lon1 = data1['lon']
time1 = data1['mtime']
pco2_1  = data1['pCO2']

pco2_1 = pco2_1.where(pco2_1 != 0) 

lat2 = data2['LAT']
lon2 = data2['LON']
time2 = data2['TIME']
pco2_2 = data2['pCO2_Original']

#%%

clim_as = pco2_0.sel(lat = slice(0,30),lon = slice(38,78)).mean(dim=('lat','lon'))

incois_as = pco2_2.sel(LAT = slice(0,30),LON = slice(38,78),TIME = slice('1980-01-01','2007-12-31'))

socat_as = pco2_1.sel(lat = slice(0,30),lon = slice(38,78),mtime = slice('1980-01-01','2007-12-31'))
#%%
socat_as = socat_as.mean(dim=('lat','lon'))
socat_mon_as = socat_as.groupby('mtime.month').mean(dim = ('mtime'))


#%%
incois_as = incois_as.mean(dim = ('LAT','LON'))
incois_mon_as = incois_as.groupby('TIME.month').mean(dim = ('TIME'))
#%%
corr = xr.corr(incois_mon_as,clim_as)
#%% AS
import matplotlib.pyplot as plt
month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']


plt.figure(figsize=(15,10), dpi=300)
plt.plot(month,clim_as, marker='o', linestyle='-', color='black',label = 'climatology')
plt.plot(month,socat_mon_as, marker='o', linestyle='-', color='green', label = ' SOCATv2024')  
plt.plot(month,incois_mon_as, marker='o', linestyle='-', color='blue',label = 'INCOIS_IBR_ROMS')  
 
plt.xlabel("Month")
plt.ylabel('\u00b5atm')
plt.title("Climatology of $pCO_2$ in Arabian Sea")
plt.grid(True)
plt.legend()
plt.show()

#%% BoB

clim_bob = pco2_0.sel(lat = slice(0,30),lon = slice(78,110)).mean(dim=('lat','lon'))

incois_bob = pco2_2.sel(LAT = slice(0,30),LON = slice(78,110),TIME = slice('1980-01-01','2007-12-31'))

socat_bob = pco2_1.sel(lat = slice(0,30),lon = slice(78,110),mtime = slice('1980-01-01','2007-12-31'))

socat_bob = socat_bob.mean(dim = ('lat','lon'))

socat_mon_bob = socat_bob.groupby('mtime.month').mean(dim = ('mtime'))

incois_bob = incois_bob.mean(dim=('LAT','LON'))

incois_mon_bob = incois_bob.groupby('TIME.month').mean(dim = ('TIME'))

#%% 

plt.figure(figsize=(15,10), dpi=300)
plt.plot(month,clim_bob, marker='o', linestyle='-', color='black',label = 'climatology')
plt.plot(month,socat_mon_bob, marker='o', linestyle='-', color='green', label = ' SOCATv2024')  
plt.plot(month,incois_mon_bob, marker='o', linestyle='-', color='blue',label = 'INCOIS_IBR_ROMS')  
 
plt.xlabel("Month")
plt.ylabel('\u00b5atm')
plt.title("Climatology of $pCO_2$ in Bay of Bengal")
plt.grid(True)
plt.legend()
plt.show()

#%% EIO

clim_eio = pco2_0.sel(lat = slice(-18,0),lon = slice(30,110)).mean(dim=('lat','lon'))

incois_eio = pco2_2.sel(LAT = slice(-18,0),LON = slice(30,110),TIME = slice('1980-01-01','2007-12-31'))

socat_eio = pco2_1.sel(lat = slice(-18,0),lon = slice(30,110),mtime = slice('1980-01-01','2007-12-31'))

socat_eio = socat_eio.mean(dim = ('lat','lon'))

socat_mon_eio = socat_eio.groupby('mtime.month').mean(dim = ('mtime'))

incois_eio= incois_eio.mean(dim=('LAT','LON'))

incois_mon_eio= incois_eio.groupby('TIME.month').mean(dim = ('TIME'))

#%% 

plt.figure(figsize=(15,10), dpi=300)
plt.plot(month,clim_eio, marker='o', linestyle='-', color='black',label = 'climatology')
plt.plot(month,socat_mon_eio, marker='o', linestyle='-', color='green', label = ' SOCATv2024')  
plt.plot(month,incois_mon_eio, marker='o', linestyle='-', color='blue',label = 'INCOIS_IBR_ROMS')  
 
plt.xlabel("Month")
plt.ylabel('\u00b5atm')
plt.title("Climatology of $pCO_2$ in Equatorial Indian Ocean")
plt.grid(True)
plt.legend()
plt.show()

#%% io

clim_io = pco2_0.sel(lat = slice(-30,30),lon = slice(30,110)).mean(dim=('lat','lon'))

incois_io = pco2_2.sel(LAT = slice(-30,30),LON = slice(30,110),TIME = slice('1980-01-01','2007-12-31'))

socat_io = pco2_1.sel(lat = slice(-30,30),lon = slice(30,110),mtime = slice('1980-01-01','2007-12-31'))

socat_io = socat_io.mean(dim = ('lat','lon'))

socat_mon_io = socat_io.groupby('mtime.month').mean(dim = ('mtime'))

incois_io= incois_io.mean(dim=('LAT','LON'))

incois_mon_io= incois_io.groupby('TIME.month').mean(dim = ('TIME'))


#%% Indian Ocean

fig, axs = plt.subplots(2, 2, figsize=(20, 15), dpi=300)  # 2x2 grid
# Individual plots
axs[0, 0].plot(month, clim_io, marker='o', linestyle='-', color='black', label='Climatology')
axs[0, 0].set_title('Indian_ocean - Climatology')
axs[0, 0].set_xlabel('Month')
axs[0, 0].set_ylabel('\u00b5atm')
axs[0, 0].grid(True)
axs[0, 0].legend()

axs[0, 1].plot(month, socat_mon_io, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[0, 1].set_title('Indian_ocean - SOCATv2024')
axs[0, 1].set_xlabel('Month')
axs[0, 1].set_ylabel('\u00b5atm')
axs[0, 1].grid(True)
axs[0, 1].legend()

axs[1, 0].plot(month, incois_mon_io, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 0].set_title('Indian_ocean - INCOIS_IBR_ROMS')
axs[1, 0].set_xlabel('Month')
axs[1, 0].set_ylabel('\u00b5atm')
axs[1, 0].grid(True)
axs[1, 0].legend()

# Combined plot
axs[1, 1].plot(month, clim_io, marker='o', linestyle='-', color='black', label='Climatology')
axs[1, 1].plot(month, socat_mon_io, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[1, 1].plot(month, incois_mon_io, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 1].set_title('Climatology of $pCO_2$ in Indian Ocean')
axs[1, 1].set_xlabel('Month')
axs[1, 1].set_ylabel('\u00b5atm')
axs[1, 1].set_ylim(200,450)
axs[1, 1].grid(True)
axs[1, 1].legend()

plt.tight_layout()
plt.show()

#%%

fig, axs = plt.subplots(2, 2, figsize=(20, 15), dpi=300)  # 2x2 grid

# --- Individual plots ---
axs[0, 0].plot(month, clim_eio, marker='o', linestyle='-', color='black', label='Climatology')
axs[0, 0].set_title('EIO - Climatology')
axs[0, 0].set_xlabel('Month')
axs[0, 0].set_ylabel('\u00b5atm')
axs[0, 0].grid(True)
axs[0, 0].legend()

axs[0, 1].plot(month, socat_mon_eio, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[0, 1].set_title('EIO - SOCATv2024')
axs[0, 1].set_xlabel('Month')
axs[0, 1].set_ylabel('\u00b5atm')
axs[0, 1].grid(True)
axs[0, 1].legend()

axs[1, 0].plot(month, incois_mon_eio, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 0].set_title('EIO - INCOIS_IBR_ROMS')
axs[1, 0].set_xlabel('Month')
axs[1, 0].set_ylabel('\u00b5atm')
axs[1, 0].grid(True)
axs[1, 0].legend()

# --- Combined plot ---
axs[1, 1].plot(month, clim_eio, marker='o', linestyle='-', color='black', label='Climatology')
axs[1, 1].plot(month, socat_mon_eio, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[1, 1].plot(month, incois_mon_eio, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 1].set_title('Climatology of $pCO_2$ in Equatorial Indian Ocean')
axs[1, 1].set_xlabel('Month')
axs[1, 1].set_ylabel('\u00b5atm')
axs[1, 1].set_ylim(200,450)
axs[1, 1].grid(True)
axs[1, 1].legend()

plt.tight_layout()
plt.show()
#%%
fig, axs = plt.subplots(2, 2, figsize=(18, 12), dpi=300)

# --- Individual Plots ---
# Climatology
axs[0, 0].plot(month, clim_as, marker='o', linestyle='-', color='black', label='Climatology')
axs[0, 0].set_title('Arabian Sea - Climatology')
axs[0, 0].set_xlabel('Month')
axs[0, 0].set_ylabel('\u00b5atm')
axs[0, 0].grid(True)
axs[0, 0].legend()

# SOCATv2024
axs[0, 1].plot(month, socat_mon_as, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[0, 1].set_title('Arabian Sea - SOCATv2024')
axs[0, 1].set_xlabel('Month')
axs[0, 1].set_ylabel('\u00b5atm')
axs[0, 1].grid(True)
axs[0, 1].legend()

# INCOIS_IBR_ROMS
axs[1, 0].plot(month, incois_mon_as, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 0].set_title('Arabian Sea - INCOIS_IBR_ROMS')
axs[1, 0].set_xlabel('Month')
axs[1, 0].set_ylabel('\u00b5atm')
axs[1, 0].grid(True)
axs[1, 0].legend()

# --- Combined Plot ---
axs[1, 1].plot(month, clim_as, marker='o', linestyle='-', color='black', label='Climatology')
axs[1, 1].plot(month, socat_mon_as, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[1, 1].plot(month, incois_mon_as, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 1].set_title('Climatology of $pCO_2$ in Arabian Sea')
axs[1, 1].set_xlabel('Month')
axs[1, 1].set_ylabel('\u00b5atm')
axs[1, 1].set_ylim(200,450)
axs[1, 1].grid(True)
axs[1, 1].legend()

plt.tight_layout()
#%%

fig, axs = plt.subplots(2, 2, figsize=(18, 12), dpi=300)


# Climatology
axs[0, 0].plot(month, clim_bob, marker='o', linestyle='-', color='black', label='Climatology')
axs[0, 0].set_title('Bay of Bengal - Climatology')
axs[0, 0].set_xlabel('Month')
axs[0, 0].set_ylabel('\u00b5atm')
axs[0, 0].grid(True)
axs[0, 0].legend()

# SOCATv2024
axs[0, 1].plot(month, socat_mon_bob, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[0, 1].set_title('Bay of Bengal - SOCATv2024')
axs[0, 1].set_xlabel('Month')
axs[0, 1].set_ylabel('\u00b5atm')
axs[0, 1].grid(True)
axs[0, 1].legend()

# INCOIS_IBR_ROMS
axs[1, 0].plot(month, incois_mon_bob, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 0].set_title('Bay of Bengal - INCOIS_IBR_ROMS')
axs[1, 0].set_xlabel('Month')
axs[1, 0].set_ylabel('\u00b5atm')
axs[1, 0].grid(True)
axs[1, 0].legend()

# --- Combined Plot ---
axs[1, 1].plot(month, clim_bob, marker='o', linestyle='-', color='black', label='Climatology')
axs[1, 1].plot(month, socat_mon_bob, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[1, 1].plot(month, incois_mon_bob, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 1].set_title('Climatology of $pCO_2$ in Bay of Bengal')
axs[1, 1].set_xlabel('Month')
axs[1, 1].set_ylabel('\u00b5atm')
axs[1, 1].set_ylim(200,450)
axs[1, 1].grid(True)
axs[1, 1].legend()
plt.tight_layout()
plt.show()
#%%
clim_sio = pco2_0.sel(lat = slice(-30,0),lon = slice(30,110)).mean(dim=('lat','lon'))

incois_sio = pco2_2.sel(LAT = slice(-30,0),LON = slice(30,110))

socat_sio = pco2_1.sel(lat = slice(-30,0),lon = slice(30,110),mtime = slice('1980-01-01','2019-12-31'))

socat_sio = socat_sio.mean(dim = ('lat','lon'))

socat_mon_sio = socat_sio.groupby('mtime.month').mean(dim = ('mtime'))

incois_sio= incois_sio.mean(dim=('LAT','LON'))

incois_mon_sio= incois_sio.groupby('TIME.month').mean(dim = ('TIME'))


#%%

fig, axs = plt.subplots(2, 2, figsize=(18, 12), dpi=300)

# --- Individual Plots ---
# Climatology
axs[0, 0].plot(month, clim_sio, marker='o', linestyle='-', color='black', label='Climatology')
axs[0, 0].set_title('Southern Indian Ocean - Climatology')
axs[0, 0].set_xlabel('Month')
axs[0, 0].set_ylabel('\u00b5atm')
axs[0, 0].grid(True)
axs[0, 0].legend()

# SOCATv2024
axs[0, 1].plot(month, socat_mon_sio, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[0, 1].set_title('Southern Indian Ocean - SOCATv2024')
axs[0, 1].set_xlabel('Month')
axs[0, 1].set_ylabel('\u00b5atm')
axs[0, 1].grid(True)
axs[0, 1].legend()

# INCOIS_IBR_ROMS
axs[1, 0].plot(month, incois_mon_sio, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 0].set_title('Southern Indian Ocean - INCOIS_IBR_ROMS')
axs[1, 0].set_xlabel('Month')
axs[1, 0].set_ylabel('\u00b5atm')
axs[1, 0].grid(True)
axs[1, 0].legend()

# --- Combined Plot ---
axs[1, 1].plot(month, clim_sio, marker='o', linestyle='-', color='black', label='Climatology')
axs[1, 1].plot(month, socat_mon_sio, marker='o', linestyle='-', color='green', label='SOCATv2024')
axs[1, 1].plot(month, incois_mon_sio, marker='o', linestyle='-', color='blue', label='INCOIS_IBR_ROMS')
axs[1, 1].set_title('Climatology of $pCO_2$ in Southern Indian Ocean')
axs[1, 1].set_xlabel('Month')
axs[1, 1].set_ylabel('\u00b5atm')
axs[1, 1].set_ylim(200,450)
axs[1, 1].grid(True)
axs[1, 1].legend()

plt.tight_layout()
plt.show()
#%%
import numpy as np

clim_as_std = np.std(clim_as)
socat_as_std = np.std(socat_mon_as)
incois_as_std = np.std(incois_mon_as)
#%%  EOF
