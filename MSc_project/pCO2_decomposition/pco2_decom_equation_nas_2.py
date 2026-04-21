#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 11:55:36 2026

@author: bobco-08
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


file = '/home/bobco-08/24cl05012/CO2/data/data_1/cmems/bio_var/cmems_obs-mob_glo_bgc-car_my_irr-i_multi-vars_35.12E-119.88E_29.88S-29.88N_1994-01-01-2024-12-01.nc'

data = xr.open_dataset(file)
data = data.rename({'latitude':'lat', 'longitude': 'lon'})
lat = data['lat']
lon =data['lon']
time = data['time']
dic = data['tco2'].sel(lat = slice(22.5,28),lon = slice(56,70))
alk = data['talk'].sel(lat = slice(22.5,28),lon = slice(56,70))
pco2 = data['spco2'].sel(lat = slice(22.5,28),lon = slice(56,70))

data2 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/cmems/phy_var/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_multi-vars_35.00E-120.00E_30.00S-30.00N_0.51m_1993-01-01-2024-12-01.nc')
data2 = data2.rename({'latitude':'lat', 'longitude': 'lon'})
sst = data2['thetao_oras'].isel(depth=0)
sss = data2['so_oras'].isel(depth=0)
sst = sst.sel(lat = slice(22.5,28),lon = slice(56,70),time = slice('1994-01-01','2024-12-01'))
sss = sss.sel(lat = slice(22.5,28),lon = slice(56,70),time = slice('1994-01-01','2024-12-01'))


#%%
pco2_reg = pco2.mean(dim=('lat', 'lon'))
pco2_annual = pco2_reg.groupby('time.year').mean(dim = 'time')

pco2_jja = pco2_reg.sel(time=pco2_reg['time'].dt.month.isin([7,8,9]))
pco2_jja_grouped = pco2_jja.groupby('time.year')
pco2_jja_y = pco2_jja_grouped.mean()

pco2_djf = pco2_reg.resample(time='QS-DEC').mean(dim='time')  
pco2_djf = pco2_djf.isel(time=slice(1,None))                 
pco2_djf = pco2_djf.sel(time=pco2_djf.time.dt.month == 12)
pco2_djf['year'] = pco2_djf['time.year'] + 1

pco2_djf_y = pco2_djf.groupby('year').mean()

delta_pco2 = (pco2_jja_y - pco2_djf_y)


#%% t-term


sst_reg = sst.mean(dim=['lat', 'lon'])

sst_jja = sst_reg.sel(time=sst_reg['time'].dt.month.isin([7,8,9]))
sst_jja_grouped = sst_jja.groupby('time.year')
sst_jja_y = sst_jja_grouped.mean()

sst_djf = sst_reg.resample(time='QS-DEC').mean(dim='time')  
sst_djf = sst_djf.isel(time=slice(1,None))                 
sst_djf = sst_djf.sel(time=sst_djf.time.dt.month == 12)
sst_djf['year'] = sst_djf['time.year'] + 1

sst_djf_y = sst_djf.groupby('year').mean()

delta_T = (sst_jja_y-sst_djf_y)


T_term = 2 * pco2_annual * (np.exp(0.0423 * (delta_T / 2)) - 1)


#%% dic_term


dic_reg = dic.mean(dim=['lat', 'lon'])
dic_annual = dic_reg.groupby('time.year').mean(dim='time')



dic_jja = dic_reg.sel(time=dic_reg['time'].dt.month.isin([7,8,9]))
dic_jja_grouped = dic_jja.groupby('time.year')
dic_jja_y = dic_jja_grouped.mean()



dic_djf = dic_reg.resample(time='QS-DEC').mean(dim='time')  
dic_djf = dic_djf.isel(time=slice(1,None))                 
dic_djf = dic_djf.sel(time=dic_djf.time.dt.month == 12)
dic_djf['year'] = dic_djf['time.year'] + 1

dic_djf_y = dic_djf.groupby('year').mean()




delta_dic = (dic_jja_y - dic_djf_y)

dic_term = (9.5 * (pco2_annual/dic_annual))* delta_dic


#%% alk_term


alk_reg = alk.mean(dim=['lat', 'lon'])
alk_annual = alk_reg.groupby('time.year').mean(dim='time')




alk_jja = alk_reg.sel(time=alk_reg['time'].dt.month.isin([7,8,9]))
alk_jja_grouped = alk_jja.groupby('time.year')
alk_jja_y = alk_jja_grouped.mean()



alk_djf = alk_reg.resample(time='QS-DEC').mean(dim='time')  
alk_djf = alk_djf.isel(time=slice(1,None))                 
alk_djf = alk_djf.sel(time=alk_djf.time.dt.month == 12)
alk_djf['year'] = alk_djf['time.year'] + 1

alk_djf_y = alk_djf.groupby('year').mean()


delta_alk = (alk_jja_y - alk_djf_y)


alk_term = (- 8.9 * (pco2_annual/alk_annual))* delta_alk


#%% sal_term

sal_reg = sss.mean(dim=['lat', 'lon'])
sal_annual = sal_reg.groupby('time.year').mean(dim='time')



sal_jja = sal_reg.sel(time=sal_reg['time'].dt.month.isin([7,8,9]))
sal_jja_grouped = sal_jja.groupby('time.year')
sal_jja_y = sal_jja_grouped.mean()



sal_djf = sal_reg.resample(time='QS-DEC').mean(dim='time')  
sal_djf = sal_djf.isel(time=slice(1,None))                 
sal_djf = sal_djf.sel(time=sal_djf.time.dt.month == 12)
sal_djf['year'] = sal_djf['time.year'] + 1

sal_djf_y = sal_djf.groupby('year').mean()

# JJAS – DJF difference
delta_sal = (sal_jja_y - sal_djf_y)


sal_term = (0.026* pco2_annual)* delta_sal

#%%
 
delta_pco2_cal = T_term + dic_term + alk_term + sal_term

#%%

x1= delta_pco2_cal.values
x2 = delta_pco2.values
sq_bias = (x1-x2)**2
rmse = np.sqrt(np.mean(sq_bias))

print(f"RMSE = {rmse:.2f} µatm")






#%% Plot observed vs calculated seasonal ΔpCO2 (1994–2024)


# Align observed and calculated ΔpCO2 in time
delta_pco2, delta_pco2_cal = xr.align(delta_pco2, delta_pco2_cal)

# Extract years and values
years = delta_pco2['year'].values
dpco2_obs_eio = delta_pco2.values
dpco2_cal_eio = delta_pco2_cal.values

plt.figure(figsize=(18,6),dpi = 400)

# Observed line
plt.plot(
    years,
    dpco2_obs_eio,
    label='Obs-reconstructed ΔpCO$_2$',
    color='k',
    linewidth=2,
    marker='o'
)

# Calculated line
plt.plot(
    years,
    dpco2_cal_eio,
    label='Calculated ΔpCO$_2$',
    color='red',
    linewidth=2,
    marker='s'
)

# mask = dpco2_cal_eio > dpco2_obs_eio
# for x, y in zip(years[mask], dpco2_cal_eio[mask]):
#     plt.vlines(x, ymin=0, ymax=y, color='red', linestyle='--', alpha=0.5)


plt.xlabel('Year',fontweight ='bold',fontsize = 16)
plt.ylabel('ΔpCO$_2$ (µatm)',fontweight ='bold',fontsize = 16)
plt.title('Obs-reconstructed vs Calculated Year-wise ΔpCO$_2$ (JAS-DJF) in NAS (1994–2024)',fontweight ='bold',fontsize= 18)
plt.xlim(1995,2024)
plt.xticks(years, rotation=45,fontsize = 12)

plt.legend(frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.4)

plt.tight_layout()
plt.show()


#%%

# ---- Align variables ----
T_term, dic_term, alk_term, sal_term = xr.align(
    T_term, dic_term, alk_term, sal_term
)

years = delta_pco2['year'].values
x = np.arange(len(years))

T   = T_term.values
DIC = dic_term.values
ALK = alk_term.values
SAL = sal_term.values

plt.figure(figsize=(20,6))

# ---- Line plots of contributions ----
plt.plot(x, T,   '-o', linewidth=2.5, label='Temperature')
plt.plot(x, DIC, '-s', linewidth=2.5, label='DIC')
plt.plot(x, ALK, '-^', linewidth=2.5, label='Alkalinity')
plt.plot(x, SAL, '-d', linewidth=2.5, label='Salinity')

# Zero reference line
plt.axhline(0, color='black', linewidth=1)

# Axis formatting
plt.xticks(x[::2], years[::2], rotation=45, fontsize=14)
plt.yticks(fontsize=14)

plt.ylabel('ΔpCO$_2$ (µatm)', fontsize=18,fontweight ='bold')
plt.title('pCO$_2$ Decomposition (JAS–DJF)\nEIO', fontsize=20)

# Legend
plt.legend(ncol=4, frameon=False, fontsize=14)

# Grid
plt.grid(axis='y', linestyle='--', alpha=0.4)

plt.tight_layout()
plt.show()
#%%
mean_dT = np.mean(delta_T)
mean_dDIC = np.mean(delta_dic)
mean_dALK = np.mean(delta_alk)
mean_dS = np.mean(delta_sal)

pco2_annual_mean = pco2_annual.mean(dim=('year'))
DIC_annual_mean = dic_annual.mean(dim=('year'))
ALK_annual_mean =alk_annual.mean(dim=('year'))

T_mean   = 2 * pco2_annual_mean * (np.exp(0.0423 * (mean_dT / 2)) - 1)
DIC_mean = (9.5 * (pco2_annual_mean/DIC_annual_mean))* mean_dDIC
ALK_mean = (- 8.9 * (pco2_annual_mean/ALK_annual_mean))* mean_dALK
SAL_mean = (0.026* pco2_annual_mean)* mean_dS


cal_pco2_mean = T_mean+DIC_mean+ALK_mean+SAL_mean 
delta_pco2_mean = np.mean(pco2_jja_y - pco2_djf_y)


sigma_SAL = np.std(sal_term, ddof=1)
sigma_T = np.std(T_term, ddof=1)
sigma_DIC = np.std(dic_term, ddof=1)
sigma_ALK = np.std(alk_term, ddof=1)
sigma_cal = np.std(delta_pco2_cal, ddof=1)
sigma_obs = np.std(delta_pco2, ddof=1)




#%%


variables = ['Temperature', 'DIC', 'Alkalinity', 'Salinity','Calculated ΔpCO$_2$','Obs-recon ΔpCO$_2$']
means = [T_mean, DIC_mean, ALK_mean, SAL_mean,cal_pco2_mean,delta_pco2_mean]

colors = ['red', 'green', 'pink', 'blue','grey','lightblue']

x = np.arange(len(variables))*0.26

fig, ax = plt.subplots(figsize=(26,14),dpi = 600)

errors = [sigma_T, sigma_DIC, sigma_ALK, sigma_SAL,sigma_cal, sigma_obs]

bars = []
for i, (var, mean, col) in enumerate(zip(variables, means, colors)):
    bar = bar = ax.bar(x[i],mean,color=col,width=0.25,label=var,yerr=errors[i],capsize=12,ecolor='black',
                       error_kw={'elinewidth':4}
)
    bars.append(bar[0])
ax.set_xlim(-0.2, x[-1] + 0.2)

ax.axhline(0, color='k', linewidth=0.8)

ax.set_ylabel('Mean ΔpCO$_2$ (µatm)',fontsize=24,fontweight = 'bold',labelpad=15)
ax.set_yticklabels(ax.get_yticks(), fontsize=24)
ax.set_title('Seasonal contrast between JAS and DJF of ΔpCO$_2$ in NAS ',fontsize=25,fontweight = 'bold')

# Remove x-axis labels
ax.set_xticks([])

ax.grid(axis='y', linestyle='--', alpha=0.2)

# ---- Add values above bars ----
# ---- Add values near the x-axis ----
y_axis_offset = max((np.array(means))) * -0.96  # small vertical offset from x-axis

for i, bar in enumerate(bars):
    val = means[i]
    text_color = 'red' if val < 0 else 'black'  # red for negative, black for positive

    ax.text(
        bar.get_x() + bar.get_width()/2,  # center of the bar
        y_axis_offset,                     # fixed height near x-axis
        f'{val:.2f} µatm',                # show mean value
        ha='center',
        va='bottom',
        fontsize=24,
        color=text_color,
        fontweight='bold'
    )

ax.set_xticks(x)
ax.set_xticklabels(variables,fontsize=24,fontweight = 'bold')

plt.tight_layout()
plt.show()
#%%

# ================================
# 📊 SUMMARY STATISTICS
# ================================

def mean_std(x):
    """Return mean and std rounded to 2 decimals"""
    return np.round(x.mean().values, 2), np.round(x.std(ddof=1).values, 2)


# ----------- ANNUAL MEAN ± STD -----------
pco2_ann_mean, pco2_ann_std = mean_std(pco2_annual)
dic_ann_mean, dic_ann_std   = mean_std(dic_annual)
alk_ann_mean, alk_ann_std   = mean_std(alk_annual)
sal_ann_mean, sal_ann_std   = mean_std(sal_annual)
sst_ann_mean, sst_ann_std   = mean_std(sst_reg.groupby('time.year').mean())


# ----------- SEASONAL AMPLITUDE (JAS–DJF) -----------
pco2_amp_mean, pco2_amp_std = mean_std(delta_pco2)
dic_amp_mean, dic_amp_std   = mean_std(delta_dic)
alk_amp_mean, alk_amp_std   = mean_std(delta_alk)
sal_amp_mean, sal_amp_std   = mean_std(delta_sal)
sst_amp_mean, sst_amp_std   = mean_std(delta_T)


# ----------- ΔpCO2 (CALCULATED & OBSERVED) -----------
obs_mean, obs_std = mean_std(delta_pco2)


# ================================
# 🖨️ PRINT RESULTS
# ================================

print("\n===== ANNUAL MEAN ± STD =====")
print(f"pCO2  : {pco2_ann_mean} ± {pco2_ann_std} µatm")
print(f"DIC   : {dic_ann_mean} ± {dic_ann_std}")
print(f"ALK   : {alk_ann_mean} ± {alk_ann_std}")
print(f"SAL   : {sal_ann_mean} ± {sal_ann_std}")
print(f"SST   : {sst_ann_mean} ± {sst_ann_std} °C")


print("\n===== SEASONAL AMPLITUDE (JAS–DJF) =====")
print(f"ΔpCO2 : {pco2_amp_mean} ± {pco2_amp_std} µatm")
print(f"ΔDIC  : {dic_amp_mean} ± {dic_amp_std}")
print(f"ΔALK  : {alk_amp_mean} ± {alk_amp_std}")
print(f"ΔSAL  : {sal_amp_mean} ± {sal_amp_std}")
print(f"ΔSST  : {sst_amp_mean} ± {sst_amp_std} °C")


print("\n===== ΔpCO2 COMPARISON =====")
print(f"Observed ΔpCO2   : {obs_mean} ± {obs_std} µatm")

# ----------- STD OF GENERATED ΔpCO2 (SEPARATE) -----------

sigma_sss = np.round(np.std(sal_term, ddof=1).values, 2)
sigma_sst   = np.round(np.std(T_term, ddof=1).values, 2)
sigma_dic = np.round(np.std(dic_term, ddof=1).values, 2)
sigma_alk = np.round(np.std(alk_term, ddof=1).values, 2)
sigma_cal_pco2 = np.round(np.std(delta_pco2_cal, ddof=1).values, 2)
sigma_obs_pco2 = np.round(np.std(delta_pco2, ddof=1).values, 2)


print("\n===== STD OF GENERATED ΔpCO2 CONTRIBUTIONS =====")
print(f"Temperature (σ) : {sigma_sst} µatm")
print(f"DIC (σ)         : {sigma_dic} µatm")
print(f"Alkalinity (σ)  : {sigma_alk} µatm")
print(f"Salinity (σ)    : {sigma_sss} µatm")

print("\n===== STD OF TOTAL ΔpCO2 =====")
print(f"Calculated (σ) : {sigma_cal_pco2} µatm")
print(f"Observed (σ)   : {sigma_obs_pco2} µatm")