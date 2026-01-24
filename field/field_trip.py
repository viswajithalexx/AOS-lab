#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 10:56:58 2026

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

file_1 = '/home/bobco-08/24cl05012/field trip/ctd_profiles_combined.csv'
file_2 = '/home/bobco-08/24cl05012/field trip/ctd_profiles_combined-1.csv'

data_1 =  pd.read_csv(file_1)

data_2 = pd.read_csv(file_2)

#%%

para_005 = data_1.loc[0:196,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)','lat','lon','cast_time_local']]


para_1_005 = para_005.loc[0:8,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
para_12_005 = para_005.loc[9:17,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

loc_1 = para_005.loc[0:8,['cast_time_local','lat','lon']]

avg_para_1_005 = (para_1_005.reset_index(drop=True)+ para_12_005.reset_index(drop=True))/2

avg_1_005 = pd.concat([loc_1,avg_para_1_005,],axis =1)

#%%
para_2_005 = para_005.loc[18:25,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
para_22_005 = para_005.loc[26:33,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

loc_2 = para_005.loc[18:25,['cast_time_local','lat','lon']].reset_index(drop =True)

   
para_23_005 = para_005.loc[26:33,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_24_005 = para_005.loc[34:41,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_25_005 = para_005.loc[34:41,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_26_005 = para_005.loc[50:57,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_27_005 = para_005.loc[58:65,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]



avg_para_2_005 = (para_2_005.reset_index(drop=True)+ para_22_005.reset_index(drop=True)+para_23_005.reset_index(drop=True)+para_24_005.reset_index(drop=True)+para_25_005.reset_index(drop=True)+para_26_005.reset_index(drop=True)+para_27_005.reset_index(drop=True))/7

avg_2_005 = pd.concat([loc_2,avg_para_2_005,],axis =1)

#%%

para_3_005 = para_005.loc[66:70,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
para_32_005 = para_005.loc[71:75,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
para_33_005 = para_005.loc[76:80,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

loc_3 = para_005.loc[66:70,['cast_time_local','lat','lon']].reset_index(drop =True)


avg_para_3_005 = (para_3_005.reset_index(drop=True)+ para_32_005.reset_index(drop=True)+para_33_005.reset_index(drop=True))/3

avg_3_005 = pd.concat([loc_3,avg_para_3_005,],axis =1)


#%%

para_4_005 = para_005.loc[81:87,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
para_42_005 = para_005.loc[88:94,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
para_43_005 = para_005.loc[95:101,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_44_005 = para_005.loc[102:108,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_45_005 = para_005.loc[109:115,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]


loc_4 = para_005.loc[81:87,['cast_time_local','lat','lon']].reset_index(drop =True)


avg_para_4_005 = (para_4_005.reset_index(drop=True)+ para_42_005.reset_index(drop=True)+para_43_005.reset_index(drop=True)+para_44_005.reset_index(drop=True)+para_45_005.reset_index(drop=True))/5

avg_4_005 = pd.concat([loc_4,avg_para_4_005,],axis =1)

#%%

para_5_005 = para_005.loc[116:123,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_52_005 = para_005.loc[124:131,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_53_005 = para_005.loc[132:139,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]


loc_5 = para_005.loc[116:123,['cast_time_local','lat','lon']].reset_index(drop =True)


avg_para_5_005 = (para_5_005.reset_index(drop=True)+ para_52_005.reset_index(drop=True)+para_53_005.reset_index(drop=True))/3

avg_5_005 = pd.concat([loc_5,avg_para_5_005,],axis =1)


#%%

para_6_005 = para_005.loc[140:146,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_62_005 = para_005.loc[147:153,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_63_005 = para_005.loc[154:160,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

loc_6 = para_005.loc[140:146,['cast_time_local','lat','lon']].reset_index(drop =True)

avg_para_6_005 = (para_6_005.reset_index(drop=True)+ para_62_005.reset_index(drop=True)+para_63_005.reset_index(drop=True))/3

avg_6_005 = pd.concat([loc_6,avg_para_6_005,],axis =1)

#%%

para_7_005 = para_005.loc[161:167,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_72_005 = para_005.loc[168:174,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_73_005 = para_005.loc[175:181,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

loc_7 = para_005.loc[161:167,['cast_time_local','lat','lon']].reset_index(drop =True)

avg_para_7_005 = (para_7_005.reset_index(drop=True)+ para_72_005.reset_index(drop=True)+para_73_005.reset_index(drop=True))/3

avg_7_005 = pd.concat([loc_7,avg_para_7_005,],axis =1)

#%%

para_8_005 = para_005.loc[182:186,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_82_005 = para_005.loc[187:191,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

para_83_005 = para_005.loc[192:196,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

loc_8 = para_005.loc[182:186,['cast_time_local','lat','lon']].reset_index(drop =True)

avg_para_8_005 = (para_8_005.reset_index(drop=True)+ para_82_005.reset_index(drop=True)+para_83_005.reset_index(drop=True))/3

avg_8_005 = pd.concat([loc_8,avg_para_8_005,],axis =1)

#%%

ctd_005_combined =  pd.concat([avg_1_005,avg_2_005,avg_3_005,avg_4_005,avg_5_005,avg_6_005,avg_7_005,avg_8_005]).reset_index(drop=True)

ctd_005_combined.to_csv('/home/bobco-08/24cl05012/field trip/ctd_005_results.csv', index=False)

#%% temperature 005
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd


TITLE_FS = 14
LABEL_FS = 13
TICK_FS  = 11


ctd_005_combined = pd.read_csv('/home/bobco-08/24cl05012/field trip/ctd_005_results.csv')


# 1. Select the first region
region_1 = ctd_005_combined.iloc[0:9]
region_2 = ctd_005_combined.iloc[9:17]
region_3 = ctd_005_combined.iloc[17:22]
region_4 = ctd_005_combined.iloc[22:29]
region_5 = ctd_005_combined.iloc[29:37]
region_6 = ctd_005_combined.iloc[37:44]
region_7 = ctd_005_combined.iloc[44:51]
region_8 = ctd_005_combined.iloc[51:56]

fig, axs = plt.subplots(2, 4, figsize=(20,15), sharey=True)

# -------- Region 1 --------
axs[0,0].plot(region_1['Temperature (Celsius)'], region_1['Depth (Meter)'], 'ro-')
axs[0,0].set_title('19.664N , 85.2138E (2025-12-09 14:50:27)', pad=18,fontsize=TITLE_FS)
axs[0,0].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[0,0].set_ylabel('Depth (m)', fontsize=LABEL_FS)
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(region_2['Temperature (Celsius)'], region_2['Depth (Meter)'], 'ro-')
axs[0,1].set_title('19.6903N , 85.2101E (2025-12-09 16:14:38)', pad=18,fontsize=TITLE_FS)
axs[0,1].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(region_3['Temperature (Celsius)'], region_3['Depth (Meter)'], 'ro-')
axs[0,2].set_title('19.7159N , 85.1956E (2025-12-09 16:56:03)', pad=18,fontsize=TITLE_FS)
axs[0,2].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(region_4['Temperature (Celsius)'], region_4['Depth (Meter)'], 'ro-')
axs[0,3].set_title('19.7266N , 85.2435E (2025-12-10 09:33:05)', pad=18,fontsize=TITLE_FS)
axs[0,3].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(region_5['Temperature (Celsius)'], region_5['Depth (Meter)'], 'ro-')
axs[1,0].set_title('19.729N , 85.2755E (2025-12-10 10:18:35)', pad=18,fontsize=TITLE_FS)
axs[1,0].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[1,0].set_ylabel('Depth (m)', fontsize=LABEL_FS)
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(region_6['Temperature (Celsius)'], region_6['Depth (Meter)'], 'ro-')
axs[1,1].set_title('19.7276N , 85.3145E (2025-12-10 11:18:26)', pad=18,fontsize=TITLE_FS)
axs[1,1].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(region_7['Temperature (Celsius)'], region_7['Depth (Meter)'], 'ro-')
axs[1,2].set_title('19.7234N , 85.2187E (2025-12-10 12:48:51)', pad=18,fontsize=TITLE_FS)
axs[1,2].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(region_8['Temperature (Celsius)'], region_8['Depth (Meter)'], 'ro-')
axs[1,3].set_title('19.7161N , 85.196E (2025-12-10 13:41:05)', pad=18,fontsize=TITLE_FS)
axs[1,3].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[1,3].grid(True, alpha=0.3)

# ---------- Common axis formatting ----------
temp_ticks = np.arange(21.5, 24.5, 0.5)
dep_ticks = np.arange(0.15, 2.65, 0.15)
for ax in axs.flat:
    ax.set_xlim(21.5, 24)
    ax.set_ylim(2.6, 0.15)
    ax.set_xticks(temp_ticks)
    ax.set_yticks(dep_ticks )
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', which='major', labelsize=TICK_FS)
plt.tight_layout()
plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_005_temp.png",dpi=500,bbox_inches="tight")
plt.show()



#%% salinity 0005

import numpy as np
import matplotlib.pyplot as plt

TITLE_FS = 14
LABEL_FS = 13
TICK_FS  = 11

fig, axs = plt.subplots(2, 4, figsize=(20,15), sharey=True)

# -------- Region 1 --------
axs[0,0].plot(region_1['Salinity (Practical Salinity Scale)'], region_1['Depth (Meter)'], 'bo-')
axs[0,0].set_title('19.664N , 85.2138E (2025-12-09 14:50:27)', pad=18, fontsize=TITLE_FS)
axs[0,0].set_xlabel('PSU', fontsize=LABEL_FS)
axs[0,0].set_ylabel('Depth (m)', fontsize=LABEL_FS)
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(region_2['Salinity (Practical Salinity Scale)'], region_2['Depth (Meter)'], 'bo-')
axs[0,1].set_title('19.6903N , 85.2101E (2025-12-09 16:14:38)', pad=18, fontsize=TITLE_FS)
axs[0,1].set_xlabel('PSU', fontsize=LABEL_FS)
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(region_3['Salinity (Practical Salinity Scale)'], region_3['Depth (Meter)'], 'bo-')
axs[0,2].set_title('19.7159N , 85.1956E (2025-12-09 16:56:03)', pad=18, fontsize=TITLE_FS)
axs[0,2].set_xlabel('PSU', fontsize=LABEL_FS)
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(region_4['Salinity (Practical Salinity Scale)'], region_4['Depth (Meter)'], 'bo-')
axs[0,3].set_title('19.7266N , 85.2435E (2025-12-10 09:33:05)', pad=18, fontsize=TITLE_FS)
axs[0,3].set_xlabel('PSU', fontsize=LABEL_FS)
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(region_5['Salinity (Practical Salinity Scale)'], region_5['Depth (Meter)'], 'bo-')
axs[1,0].set_title('19.729N , 85.2755E (2025-12-10 10:18:35)', pad=18, fontsize=TITLE_FS)
axs[1,0].set_xlabel('PSU', fontsize=LABEL_FS)
axs[1,0].set_ylabel('Depth (m)', fontsize=LABEL_FS)
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(region_6['Salinity (Practical Salinity Scale)'], region_6['Depth (Meter)'], 'bo-')
axs[1,1].set_title('19.7276N , 85.3145E (2025-12-10 11:18:26)', pad=18, fontsize=TITLE_FS)
axs[1,1].set_xlabel('PSU', fontsize=LABEL_FS)
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(region_7['Salinity (Practical Salinity Scale)'], region_7['Depth (Meter)'], 'bo-')
axs[1,2].set_title('19.7234N , 85.2187E (2025-12-10 12:48:51)', pad=18, fontsize=TITLE_FS)
axs[1,2].set_xlabel('PSU', fontsize=LABEL_FS)
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(region_8['Salinity (Practical Salinity Scale)'], region_8['Depth (Meter)'], 'bo-')
axs[1,3].set_title('19.7161N , 85.196E (2025-12-10 13:41:05)', pad=18, fontsize=TITLE_FS)
axs[1,3].set_xlabel('PSU', fontsize=LABEL_FS)
axs[1,3].grid(True, alpha=0.3)

# ---------- Common axis formatting ----------
sal_ticks = np.arange(7.6, 9.6, 0.2)
dep_ticks = np.arange(0.15, 2.65, 0.15)

for ax in axs.flat:
    ax.set_xlim(7.6, 9.6)
    ax.set_ylim(2.6, 0.15)
    ax.set_xticks(sal_ticks)
    ax.set_yticks(dep_ticks)
    ax.tick_params(axis='both', labelsize=TICK_FS)
    ax.minorticks_on()
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

plt.tight_layout()
plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_005_psu.png",
            dpi=500, bbox_inches="tight")
plt.show()

#%% sound 005

import numpy as np
import matplotlib.pyplot as plt

TITLE_FS = 14
LABEL_FS = 13
TICK_FS  = 11

fig, axs = plt.subplots(2, 4, figsize=(20,15), sharey=True)

# -------- Region 1 --------
axs[0,0].plot(region_1['Sound velocity (Meters per Second)'], region_1['Depth (Meter)'], 'o-')
axs[0,0].set_title('19.664N , 85.2138E (2025-12-09 14:50:27)', pad=18, fontsize=TITLE_FS)
axs[0,0].set_xlabel('Sound velocity (m/s)', fontsize=LABEL_FS)
axs[0,0].set_ylabel('Depth (m)', fontsize=LABEL_FS)
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(region_2['Sound velocity (Meters per Second)'], region_2['Depth (Meter)'], 'o-')
axs[0,1].set_title('19.6903N , 85.2101E (2025-12-09 16:14:38)', pad=18, fontsize=TITLE_FS)
axs[0,1].set_xlabel('Sound velocity (m/s)', fontsize=LABEL_FS)
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(region_3['Sound velocity (Meters per Second)'], region_3['Depth (Meter)'], 'o-')
axs[0,2].set_title('19.7159N , 85.1956E (2025-12-09 16:56:03)', pad=18, fontsize=TITLE_FS)
axs[0,2].set_xlabel('Sound velocity (m/s)', fontsize=LABEL_FS)
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(region_4['Sound velocity (Meters per Second)'], region_4['Depth (Meter)'], 'o-')
axs[0,3].set_title('19.7266N , 85.2435E (2025-12-10 09:33:05)', pad=18, fontsize=TITLE_FS)
axs[0,3].set_xlabel('Sound velocity (m/s)', fontsize=LABEL_FS)
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(region_5['Sound velocity (Meters per Second)'], region_5['Depth (Meter)'], 'o-')
axs[1,0].set_title('19.729N , 85.2755E (2025-12-10 10:18:35)', pad=18, fontsize=TITLE_FS)
axs[1,0].set_xlabel('Sound velocity (m/s)', fontsize=LABEL_FS)
axs[1,0].set_ylabel('Depth (m)', fontsize=LABEL_FS)
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(region_6['Sound velocity (Meters per Second)'], region_6['Depth (Meter)'], 'o-')
axs[1,1].set_title('19.7276N , 85.3145E (2025-12-10 11:18:26)', pad=18, fontsize=TITLE_FS)
axs[1,1].set_xlabel('Sound velocity (m/s)', fontsize=LABEL_FS)
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(region_7['Sound velocity (Meters per Second)'], region_7['Depth (Meter)'], 'o-')
axs[1,2].set_title('19.7234N , 85.2187E (2025-12-10 12:48:51)', pad=18, fontsize=TITLE_FS)
axs[1,2].set_xlabel('Sound velocity (m/s)', fontsize=LABEL_FS)
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(region_8['Sound velocity (Meters per Second)'], region_8['Depth (Meter)'], 'o-')
axs[1,3].set_title('19.7161N , 85.196E (2025-12-10 13:41:05)', pad=18, fontsize=TITLE_FS)
axs[1,3].set_xlabel('Sound velocity (m/s)', fontsize=LABEL_FS)
axs[1,3].grid(True, alpha=0.3)

# ---------- Common axis formatting ----------
sound_ticks = np.linspace(1496, 1504, 9)
dep_ticks = np.arange(0.15, 2.65, 0.15)

for ax in axs.flat:
    ax.set_xlim(1496, 1504)
    ax.set_ylim(2.6, 0.15)
    ax.set_xticks(sound_ticks)
    ax.set_yticks(dep_ticks)
    ax.tick_params(axis='both', labelsize=TICK_FS)
    ax.minorticks_on()
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

plt.tight_layout()
plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_005_sound_velocity.png",
            dpi=500, bbox_inches="tight")
plt.show()


#%% den 005
import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 4, figsize=(20,15), sharey=True)

# -------- Region 1 --------
axs[0,0].plot(region_1['Density (Kilograms per Cubic Meter)'], region_1['Depth (Meter)'], 'go-')
axs[0,0].set_title('19.664N , 85.2138E (2025-12-09 14:50:27)', pad=18,fontsize=TITLE_FS)
axs[0,0].set_xlabel('Density (kg/m³)',fontsize=LABEL_FS)
axs[0,0].set_ylabel('Depth (m)',fontsize=LABEL_FS)
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(region_2['Density (Kilograms per Cubic Meter)'], region_2['Depth (Meter)'], 'go-')
axs[0,1].set_title('19.6903N , 85.2101E (2025-12-09 16:14:38)', pad=18,fontsize=TITLE_FS)
axs[0,1].set_xlabel('Density (kg/m³)',fontsize=LABEL_FS)
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(region_3['Density (Kilograms per Cubic Meter)'], region_3['Depth (Meter)'], 'go-')
axs[0,2].set_title('19.7159N , 85.1956E (2025-12-09 16:56:03)', pad=18,fontsize=TITLE_FS)
axs[0,2].set_xlabel('Density (kg/m³)',fontsize=LABEL_FS)
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(region_4['Density (Kilograms per Cubic Meter)'], region_4['Depth (Meter)'], 'go-')
axs[0,3].set_title('19.7266N , 85.2435E (2025-12-10 09:33:05)', pad=18,fontsize=TITLE_FS)
axs[0,3].set_xlabel('Density (kg/m³)',fontsize=LABEL_FS)
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(region_5['Density (Kilograms per Cubic Meter)'], region_5['Depth (Meter)'], 'go-')
axs[1,0].set_title('19.729N , 85.2755E (2025-12-10 10:18:35)', pad=18,fontsize=TITLE_FS)
axs[1,0].set_xlabel('Density (kg/m³)',fontsize=LABEL_FS)
axs[1,0].set_ylabel('Depth (m)',fontsize=LABEL_FS)
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(region_6['Density (Kilograms per Cubic Meter)'], region_6['Depth (Meter)'], 'go-')
axs[1,1].set_title('19.7276N , 85.3145E (2025-12-10 11:18:26)', pad=18,fontsize=TITLE_FS)
axs[1,1].set_xlabel('Density (kg/m³)',fontsize=LABEL_FS)
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(region_7['Density (Kilograms per Cubic Meter)'], region_7['Depth (Meter)'], 'go-')
axs[1,2].set_title('19.7234N , 85.2187E (2025-12-10 12:48:51)', pad=18,fontsize=TITLE_FS)
axs[1,2].set_xlabel('Density (kg/m³)',fontsize=LABEL_FS)
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(region_8['Density (Kilograms per Cubic Meter)'], region_8['Depth (Meter)'], 'go-')
axs[1,3].set_title('19.7161N , 85.196E (2025-12-10 13:41:05)', pad=18,fontsize=TITLE_FS)
axs[1,3].set_xlabel('Density (kg/m³)',fontsize=LABEL_FS)
axs[1,3].grid(True, alpha=0.3)

#---------- Common axis formatting ----------
# You can adjust x-limits according to your density range
density_ticks = np.arange(1003,1005.5,0.5)  # example typical seawater density
dep_ticks = np.arange(0.15, 2.65, 0.15)

for ax in axs.flat:
    ax.set_xlim(1003,1005)
    ax.set_ylim(2.6, 0.15)  # invert depth axis
    ax.set_xticks(density_ticks)
    ax.set_yticks(dep_ticks)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

plt.tight_layout()
plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_005_dens.png",dpi=500,bbox_inches="tight")
plt.show()


#%% tds 005

import matplotlib.pyplot as plt
import numpy as np

# Common ticks
temp_ticks = np.arange(21.5, 24.5, 0.5)
sal_ticks = np.linspace(7.6, 9.6, 6)
dens_ticks = np.arange(1003,1005.8,0.8)
dep_ticks = np.arange(0.15, 2.75, 0.15)

regions = [region_1, region_2, region_3, region_4,
           region_5, region_6, region_7, region_8]

titles = [
    '19.664N , 85.2138E (2025-12-09 14:50:27)',
    '19.6903N , 85.2101E (2025-12-09 16:14:38)',
    '19.7159N , 85.1956E (2025-12-09 16:56:03)',
    '19.7266N , 85.2435E (2025-12-10 09:33:05)',
    '19.729N , 85.2755E (2025-12-10 10:18:35)',
    '19.7276N , 85.3145E (2025-12-10 11:18:26)',
    '19.7234N , 85.2187E (2025-12-10 12:48:51)',
    '19.7161N , 85.196E (2025-12-10 13:41:05)'
]

fig, axs = plt.subplots(2, 4, figsize=(20, 15), sharey=True)

for i, (region, ax) in enumerate(zip(regions, axs.flat)):
    # 1. PRIMARY X-AXIS: Salinity (Bottom)
    ax.plot(region['Salinity (Practical Salinity Scale)'], region['Depth (Meter)'],
            'bo-', label='Salinity')
    ax.set_xlabel('Salinity (PSU)', color='blue', fontsize=LABEL_FS)
    ax.set_xlim(sal_ticks[0], sal_ticks[-1])
    ax.set_ylim(dep_ticks[-1], dep_ticks[0])  # Invert depth
    ax.tick_params(axis='x', labelcolor='blue', labelsize=TICK_FS)

    # 2. SECOND X-AXIS: Temperature (Top)
    ax_temp = ax.twiny()
    ax_temp.plot(region['Temperature (Celsius)'], region['Depth (Meter)'],
                 'ro-', label='Temperature')
    ax_temp.set_xlabel('Temperature (°C)', color='red', fontsize=LABEL_FS)
    ax_temp.set_xlim(temp_ticks[0], temp_ticks[-1])
    ax_temp.tick_params(axis='x', labelcolor='red', labelsize=TICK_FS)

    # 3. THIRD X-AXIS: Density (Offset Top)
    ax_dens = ax.twiny()
    ax_dens.spines['top'].set_position(('outward', 40)) 
    ax_dens.plot(region['Density (Kilograms per Cubic Meter)'], region['Depth (Meter)'],
                 'go-', label='Density')
    ax_dens.set_xlabel('Density (kg/m³)', color='green', fontsize=LABEL_FS)
    ax_dens.set_xlim(dens_ticks[0], dens_ticks[-1])
    ax_dens.tick_params(axis='x', labelcolor='green', labelsize=TICK_FS)

    # General Formatting
    ax.set_ylabel('Depth (Meter)', fontsize=LABEL_FS)
    ax.tick_params(axis='y', labelsize=TICK_FS)
    ax.grid(True, alpha=0.3)
    ax.set_title(titles[i], pad=10, fontsize=TITLE_FS)  # Title font size

plt.subplots_adjust(
    left=0.05, 
    right=0.95, 
    top=0.85, 
    bottom=0.1, 
    hspace=0.7, 
    wspace=0.3
)

plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_005_tsd.png",
            dpi=500, bbox_inches="tight")
plt.show()
#%% temperature 006

ctd_006_combined = data_2.loc[0:43,['cast_time_local','lat','lon','Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]


# 1. Select the first region
reg006_1 = ctd_006_combined.iloc[0:8]
reg006_2 = ctd_006_combined.iloc[8:14]
reg006_3 = ctd_006_combined.iloc[14:17]
reg006_4 = ctd_006_combined.iloc[17:23]
reg006_5 = ctd_006_combined.iloc[23:29]
reg006_6 = ctd_006_combined.iloc[29:34]
reg006_7 = ctd_006_combined.iloc[34:40]
reg006_8 = ctd_006_combined.iloc[40:44]

fig, axs = plt.subplots(2, 4, figsize=(20,15), sharey=True)

# Font sizes
TITLE_FS = 14
LABEL_FS = 13
TICK_FS = 11

# -------- Region 1 --------
axs[0,0].plot(reg006_1['Temperature (Celsius)'], reg006_1['Depth (Meter)'], 'ro-')
axs[0,0].set_title('19.6641N , 85.2148E (2025-12-09 14:44:52)', pad=18, fontsize=TITLE_FS)
axs[0,0].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[0,0].set_ylabel('Depth (m)', fontsize=LABEL_FS)
axs[0,0].tick_params(axis='both', labelsize=TICK_FS)
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(reg006_2['Temperature (Celsius)'], reg006_2['Depth (Meter)'], 'ro-')
axs[0,1].set_title('19.6898N , 85.2097E (2025-12-09 16:15:36)', pad=18, fontsize=TITLE_FS)
axs[0,1].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[0,1].tick_params(axis='both', labelsize=TICK_FS)
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(reg006_3['Temperature (Celsius)'], reg006_3['Depth (Meter)'], 'ro-')
axs[0,2].set_title('19.7159N , 85.1958E (2025-12-09 16:54:37)', pad=18, fontsize=TITLE_FS)
axs[0,2].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].tick_params(axis='both', labelsize=TICK_FS)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(reg006_4['Temperature (Celsius)'], reg006_4['Depth (Meter)'], 'ro-')
axs[0,3].set_title('19.7271N , 85.2438E (2025-12-10 09:32:12)', pad=18, fontsize=TITLE_FS)
axs[0,3].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[0,3].tick_params(axis='both', labelsize=TICK_FS)
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(reg006_5['Temperature (Celsius)'], reg006_5['Depth (Meter)'], 'ro-')
axs[1,0].set_title('19.7297N , 85.2763E (2025-12-10 10:18:39)', pad=18, fontsize=TITLE_FS)
axs[1,0].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[1,0].set_ylabel('Depth (m)', fontsize=LABEL_FS)
axs[1,0].tick_params(axis='both', labelsize=TICK_FS)
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(reg006_6['Temperature (Celsius)'], reg006_6['Depth (Meter)'], 'ro-')
axs[1,1].set_title('19.7259N , 85.3148E (2025-12-10 11:19:17)', pad=18, fontsize=TITLE_FS)
axs[1,1].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[1,1].tick_params(axis='both', labelsize=TICK_FS)
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(reg006_7['Temperature (Celsius)'], reg006_7['Depth (Meter)'], 'ro-')
axs[1,2].set_title('19.7229N , 85.2174E (2025-12-10 12:48:02)', pad=18, fontsize=TITLE_FS)
axs[1,2].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[1,2].tick_params(axis='both', labelsize=TICK_FS)
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(reg006_8['Temperature (Celsius)'], reg006_8['Depth (Meter)'], 'ro-')
axs[1,3].set_title('19.7164N , 85.1964E (2025-12-10 13:41:04)', pad=18, fontsize=TITLE_FS)
axs[1,3].set_xlabel('Temp (°C)', fontsize=LABEL_FS)
axs[1,3].tick_params(axis='both', labelsize=TICK_FS)
axs[1,3].grid(True, alpha=0.3)

# ---------- Common axis formatting ----------
temp_ticks = np.arange(21.5, 24.5, 0.5)
dep_ticks = np.arange(0.15, 2.65, 0.15)
for ax in axs.flat:
    ax.set_xlim(21.6, 24)
    ax.set_ylim(2.6, 0.15)
    ax.set_xticks(temp_ticks)
    ax.set_yticks(dep_ticks)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

plt.tight_layout()
plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_006_temp.png",
            dpi=500, bbox_inches="tight")
plt.show()



#%% salinity 006
fig, axs = plt.subplots(2, 4, figsize=(20,15), sharey=True)
# -------- Region 1 --------
axs[0,0].plot(reg006_1['Salinity (Practical Salinity Scale)'], reg006_1['Depth (Meter)'], 'bo-')
axs[0,0].set_title('19.6641N , 85.2148E (2025-12-09 14:44:52)', pad=18,fontsize=TITLE_FS)
axs[0,0].set_xlabel('PSU',fontsize=LABEL_FS)
axs[0,0].set_ylabel('Depth (m)',fontsize=LABEL_FS)
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(reg006_2['Salinity (Practical Salinity Scale)'], reg006_2['Depth (Meter)'], 'bo-')
axs[0,1].set_title('19.6898N , 85.2097E (2025-12-09 16:15:36)', pad=18,fontsize=TITLE_FS)
axs[0,1].set_xlabel('PSU',fontsize=LABEL_FS)
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(reg006_3['Salinity (Practical Salinity Scale)'], reg006_3['Depth (Meter)'], 'bo-')
axs[0,2].set_title('19.7159N , 85.1958E (2025-12-09 16:54:37)', pad=18,fontsize=TITLE_FS)
axs[0,2].set_xlabel('PSU',fontsize=LABEL_FS)
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(reg006_4['Salinity (Practical Salinity Scale)'], reg006_4['Depth (Meter)'], 'bo-')
axs[0,3].set_title('19.7271N , 85.2438E (2025-12-10 09:32:12)', pad=18,fontsize=TITLE_FS)
axs[0,3].set_xlabel('PSU',fontsize=LABEL_FS)
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(reg006_5['Salinity (Practical Salinity Scale)'], reg006_5['Depth (Meter)'], 'bo-')
axs[1,0].set_title('19.7297N , 85.2763E (2025-12-10 10:18:39)', pad=18,fontsize=TITLE_FS)
axs[1,0].set_xlabel('PSU)',fontsize=LABEL_FS)
axs[1,0].set_ylabel('Depth (m)',fontsize=LABEL_FS)
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(reg006_6['Salinity (Practical Salinity Scale)'], reg006_6['Depth (Meter)'], 'bo-')
axs[1,1].set_title('19.7259N , 85.3148E (2025-12-10 11:19:17)', pad=18,fontsize=TITLE_FS)
axs[1,1].set_xlabel('PSU',fontsize=LABEL_FS)
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(reg006_7['Salinity (Practical Salinity Scale)'], reg006_7['Depth (Meter)'], 'bo-')
axs[1,2].set_title('19.7229N , 85.2174E (2025-12-10 12:48:02)', pad=18,fontsize=TITLE_FS)
axs[1,2].set_xlabel('PSU',fontsize=LABEL_FS)
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(reg006_8['Salinity (Practical Salinity Scale)'], reg006_8['Depth (Meter)'], 'bo-')
axs[1,3].set_title('19.7164N , 85.1964E (2025-12-10 13:41:04)', pad=18,fontsize=TITLE_FS)
axs[1,3].set_xlabel('PSU',fontsize=LABEL_FS)
axs[1,3].grid(True, alpha=0.3)

# ---------- Common axis formatting ----------
sal_ticks = np.arange(7.1, 8.05, 0.1)
dep_ticks = np.arange(0.15, 2.65, 0.15)
for ax in axs.flat:
    ax.set_xlim(7.1, 8.05)
    ax.set_ylim(2.6, 0.15)
    ax.set_xticks(sal_ticks)
    ax.set_yticks(dep_ticks )
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
plt.tight_layout()
plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_006_psu.png",dpi=500,bbox_inches="tight")
plt.show()

#%% sound  006



fig, axs = plt.subplots(2, 4, figsize=(20,15), dpi =150, sharey=True)

# -------- Region 1 --------
axs[0,0].plot(reg006_1['Sound velocity (Meters per Second)'], reg006_1['Depth (Meter)'], 'o-')
axs[0,0].set_title('19.6641N , 85.2148E (2025-12-09 14:44:52)', pad=18,fontsize=TITLE_FS)
axs[0,0].set_xlabel('Sound velocity (m/s)',fontsize=LABEL_FS)
axs[0,0].set_ylabel('Depth (m)',fontsize=LABEL_FS)
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(reg006_2['Sound velocity (Meters per Second)'], reg006_2['Depth (Meter)'], 'o-')
axs[0,1].set_title('19.6898N , 85.2097E (2025-12-09 16:15:36)', pad=18,fontsize=TITLE_FS)
axs[0,1].set_xlabel('Sound velocity (m/s)',fontsize=LABEL_FS)
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(reg006_3['Sound velocity (Meters per Second)'], reg006_3['Depth (Meter)'], 'o-')
axs[0,2].set_title('19.7159N , 85.1958E (2025-12-09 16:54:37)', pad=18,fontsize=TITLE_FS)
axs[0,2].set_xlabel('Sound velocity (m/s)',fontsize=LABEL_FS)
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(reg006_4['Sound velocity (Meters per Second)'], reg006_4['Depth (Meter)'], 'o-')
axs[0,3].set_title('19.7271N , 85.2438E (2025-12-10 09:32:12)', pad=18,fontsize=TITLE_FS)
axs[0,3].set_xlabel('Sound velocity (m/s)',fontsize=LABEL_FS)
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(reg006_5['Sound velocity (Meters per Second)'], reg006_5['Depth (Meter)'], 'o-')
axs[1,0].set_title('19.7297N , 85.2763E (2025-12-10 10:18:39)', pad=18,fontsize=TITLE_FS)
axs[1,0].set_xlabel('Sound velocity (m/s)',fontsize=LABEL_FS)
axs[1,0].set_ylabel('Depth (m)',fontsize=LABEL_FS)
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(reg006_6['Sound velocity (Meters per Second)'], reg006_6['Depth (Meter)'], 'o-')
axs[1,1].set_title('19.7259N , 85.3148E (2025-12-10 11:19:17)', pad=18,fontsize=TITLE_FS)
axs[1,1].set_xlabel('Sound velocity (m/s)',fontsize=LABEL_FS)
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(reg006_7['Sound velocity (Meters per Second)'], reg006_7['Depth (Meter)'], 'o-')
axs[1,2].set_title('19.7229N , 85.2174E (2025-12-10 12:48:02)', pad=18,fontsize=TITLE_FS)
axs[1,2].set_xlabel('Sound velocity (m/s)',fontsize=LABEL_FS)
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(reg006_8['Sound velocity (Meters per Second)'], reg006_8['Depth (Meter)'], 'o-')
axs[1,3].set_title('19.7164N , 85.1964E (2025-12-10 13:41:04)', pad=18,fontsize=TITLE_FS)
axs[1,3].set_xlabel('Sound velocity (m/s)',fontsize=LABEL_FS)
axs[1,3].grid(True, alpha=0.3)

# ---------- Common axis formatting ----------
sound_ticks = np.linspace(1494.8,1501.7, 5)
dep_ticks = np.arange(0.15, 2.65, 0.15)
for ax in axs.flat:
    ax.set_xlim(1495.4,1501.65)
    ax.set_ylim(2.6, 0.15)
    ax.set_xticks(sound_ticks)
    ax.set_yticks(dep_ticks )
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
plt.tight_layout()
plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_006_sound.png",dpi=500,bbox_inches="tight")
plt.show()

#%% density 006


fig, axs = plt.subplots(2, 4, figsize=(20,15), dpi =150, sharey=True)

# -------- Region 1 --------
axs[0,0].plot(reg006_1['Density (Kilograms per Cubic Meter)'], reg006_1['Depth (Meter)'], 'go-')
axs[0,0].set_title('19.6641N , 85.2148E (2025-12-09 14:44:52)', pad=18,fontsize=TITLE_FS)
axs[0,0].set_xlabel('Density (kg/m³)')
axs[0,0].set_ylabel('Depth (m)')

axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(reg006_2['Density (Kilograms per Cubic Meter)'], reg006_2['Depth (Meter)'], 'go-')
axs[0,1].set_title('19.6898N , 85.2097E (2025-12-09 16:15:36)', pad=18,fontsize=TITLE_FS)
axs[0,1].set_xlabel('Density (kg/m³)')
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(reg006_3['Density (Kilograms per Cubic Meter)'], reg006_3['Depth (Meter)'], 'go-')
axs[0,2].set_title('19.7159N , 85.1958E (2025-12-09 16:54:37)', pad=18,fontsize=TITLE_FS)
axs[0,2].set_xlabel('Density (kg/m³)')

axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(reg006_4['Density (Kilograms per Cubic Meter)'], reg006_4['Depth (Meter)'], 'go-')
axs[0,3].set_title('19.7271N , 85.2438E (2025-12-10 09:32:12)', pad=18,fontsize=TITLE_FS)
axs[0,3].set_xlabel('Density (kg/m³)')

axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(reg006_5['Density (Kilograms per Cubic Meter)'], reg006_5['Depth (Meter)'], 'go-')
axs[1,0].set_title('19.7297N , 85.2763E (2025-12-10 10:18:39)', pad=18,fontsize=TITLE_FS)
axs[1,0].set_xlabel('Density (kg/m³)')
axs[1,0].set_ylabel('Depth (m)')

axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(reg006_6['Density (Kilograms per Cubic Meter)'], reg006_6['Depth (Meter)'], 'go-')
axs[1,1].set_title('19.7259N , 85.3148E (2025-12-10 11:19:17)', pad=18,fontsize=TITLE_FS)
axs[1,1].set_xlabel('Density (kg/m³)')

axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(reg006_7['Density (Kilograms per Cubic Meter)'], reg006_7['Depth (Meter)'], 'go-')
axs[1,2].set_title('19.7229N , 85.2174E (2025-12-10 12:48:02)', pad=18,fontsize=TITLE_FS)
axs[1,2].set_xlabel('Density (kg/m³)')

axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(reg006_8['Density (Kilograms per Cubic Meter)'], reg006_8['Depth (Meter)'], 'go-')
axs[1,3].set_title('19.7164N , 85.1964E (2025-12-10 13:41:04)', pad=18,fontsize=TITLE_FS)
axs[1,3].set_xlabel('Density (kg/m³)')

axs[1,3].grid(True, alpha=0.3)

# ---------- Common axis formatting ----------
dens_ticks = np.arange(1003.10, 1003.56, 0.05)
dep_ticks  = np.arange(0.15, 2.65, 0.15)

for ax in axs.flat:
    ax.set_xlim(1003.10, 1003.55)
    ax.set_ylim(2.6, 0.15)   # inverted depth
    ax.set_xticks(dens_ticks)
    ax.set_yticks(dep_ticks)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

plt.tight_layout()
plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_006_den.png", dpi=500, bbox_inches="tight")
plt.show()

#%%  tds_006


# Common ticks
temp_ticks = np.linspace(21.55, 23.8, 9)
sal_ticks = np.linspace(7.15, 7.95, 8)
dens_ticks =  np.arange(1003.10, 1003.56, 0.05)
dep_ticks = np.arange(0.15, 2.75,0.15)

regions1 = [reg006_1, reg006_2, reg006_3, reg006_4,
           reg006_5, reg006_6, reg006_7, reg006_8]

titles = [
    '19.6641N , 85.2148E (2025-12-09 14:44:52)',
    '19.6898N , 85.2097E (2025-12-09 16:15:36)',
    '19.7159N , 85.1958E (2025-12-09 16:54:37)',
    '19.7271N , 85.2438E (2025-12-10 09:32:12)',
    '19.7297N , 85.2763E (2025-12-10 10:18:39)',
    '19.7259N , 85.3148E (2025-12-10 11:19:17)',
    '19.7234N , 85.2187E (2025-12-10 12:48:51)',
    '19.7161N , 85.196E (2025-12-10 13:41:05)'
]

# -------- FONT SIZES --------
TITLE_FS = 14
LABEL_FS = 13
TICK_FS  = 11

fig, axs = plt.subplots(2, 4, figsize=(20, 15), sharey=True)

for i, (region, ax) in enumerate(zip(regions1, axs.flat)):

    # ---- Salinity (bottom axis) ----
    ax.plot(region['Salinity (Practical Salinity Scale)'],
            region['Depth (Meter)'], 'bo-')
    ax.set_xlabel('Salinity (PSU)', color='blue', fontsize=LABEL_FS)
    ax.set_xlim(sal_ticks[0], sal_ticks[-1])
    ax.set_ylim(dep_ticks[-1], dep_ticks[0])
    ax.tick_params(axis='x', labelcolor='blue', labelsize=TICK_FS)
    ax.tick_params(axis='y', labelsize=TICK_FS)

    # ---- Temperature (top axis) ----
    ax_temp = ax.twiny()
    ax_temp.plot(region['Temperature (Celsius)'],
                 region['Depth (Meter)'], 'ro-')
    ax_temp.set_xlabel('Temperature (°C)', color='red', fontsize=LABEL_FS)
    ax_temp.set_xlim(temp_ticks[0], temp_ticks[-1])
    ax_temp.tick_params(axis='x', labelcolor='red', labelsize=TICK_FS)

    # ---- Density (offset top axis) ----
    ax_dens = ax.twiny()
    ax_dens.spines['top'].set_position(('outward', 40))
    ax_dens.plot(region['Density (Kilograms per Cubic Meter)'],
                 region['Depth (Meter)'], 'go-')
    ax_dens.set_xlabel('Density (kg/m³)', color='green', fontsize=LABEL_FS)
    ax_dens.set_xlim(dens_ticks[0], dens_ticks[-1])
    ax_dens.tick_params(axis='x', labelcolor='green', labelsize=TICK_FS)

    # ---- General formatting ----
    ax.set_ylabel('Depth (Meter)', fontsize=LABEL_FS)
    ax.set_title(titles[i], pad=12, fontsize=TITLE_FS)
    ax.grid(True, alpha=0.3)

plt.subplots_adjust(
    left=0.05, right=0.95, top=0.85, bottom=0.1,
    hspace=0.7, wspace=0.3
)

plt.savefig("/home/bobco-08/24cl05012/field trip/plots/ctd_006_tsd.png",
            dpi=500, bbox_inches="tight")
plt.show()


#%%
file_3 = '/home/bobco-08/24cl05012/field trip/ctd_field_seabird.xlsx'

data_3 =  pd.read_excel(file_3)

ctdsb_1 = data_3.loc[0:7,['DateTime (UTC+05:30)','Temperature (Celsius)','Conductivity (uS/cm)','Pressure (Decibar)','Oxygen (mg/L)','Chlorophyll (ug/l)','Turbidity (NTU)','Salinity (psu)','Spec Conductivity (uS/cm)','Oxygen Sat (%)']]

#%%

track_1 = ctd_005_combined[["lat", "lon"]].drop_duplicates().reset_index(drop=True)


track_2 = ctd_006_combined[["lat", "lon"]].drop_duplicates().reset_index(drop=True)

import numpy as np

# combine both tracks to get common limits
all_lon = np.concatenate([track_1["lon"], track_2["lon"]])
all_lat = np.concatenate([track_1["lat"], track_2["lat"]])

pad_lon = 0.01   # adjust if needed
pad_lat = 0.01

plt.figure()
plt.plot(track_1["lon"], track_1["lat"], "-o", label="Boat 1")
plt.plot(track_2["lon"], track_2["lat"], "-o", label="Boat 2")

plt.xlim(all_lon.min() - pad_lon, all_lon.max() + pad_lon)
plt.ylim(all_lat.min() - pad_lat, all_lat.max() + pad_lat)

plt.xlabel("Longitude (°E)")
plt.ylabel("Latitude (°N)")
plt.title("Parallel CTD Boat Tracks (Zoomed)")
plt.legend()
plt.grid(True)
plt.show()








#%% gradient (t,S,D)


t_11 = ctd_005_combined.loc[0:7,['Temperature (Celsius)']].to_numpy()
t_12 = ctd_006_combined.loc[0:7,['Temperature (Celsius)']].to_numpy()

dT1 = t_11 - t_12

from haversine import haversine
import numpy as np

D1 = np.array([
    haversine((lat1, lon1), (lat2, lon2))
    for lat1, lon1, lat2, lon2 in zip(
        ctd_005_combined.lat.iloc[:8],
        ctd_005_combined.lon.iloc[:8],
        ctd_006_combined.lat.iloc[:8],
        ctd_006_combined.lon.iloc[:8],
    )
])

grad_t1 = dT1/D1

depth1 = np.linspace(0.15, 2.14,8)

lon1 = np.linspace(85.2138,85.2148,8)

dist = np.linspace(0,0.106503,8)*1000

fig, ax1 = plt.subplots(figsize=(10,5))

pcm = ax1.contourf(lon1, depth1, grad_t1,cmap = 'viridis')
contours = ax1.contour(lon1, depth1, grad_t1, colors='black', linewidths=1)
ax1.clabel(contours, inline=True, fontsize=8, fmt="%.2f")  # Add labels to contours
ax1.invert_yaxis()
ax1.set_xlabel("Longitude (°E)")
ax1.ticklabel_format(style='plain', axis='x', useOffset=False)
ax1.set_ylabel("Depth (m)")
plt.colorbar(pcm, ax=ax1, label="dT/dy (°C/m)")

ax2 = ax1.twiny()
ax2.set_xlim(depth1.min(),depth1.max())
ax2.set_xlabel("Distance (m)")

plt.title("Horizontal Temperature Gradient Section 1")
plt.show()

#%%


t_21 = ctd_005_combined.loc[9:13,['Temperature (Celsius)']].to_numpy()
t_22 = ctd_006_combined.loc[8:12,['Temperature (Celsius)']].to_numpy()

dT2 = t_21 - t_22

from haversine import haversine
import numpy as np

D2 = np.array([
    haversine((lat1, lon1), (lat2, lon2))
    for lat1, lon1, lat2, lon2 in zip(
        ctd_005_combined.lat.iloc[9:13],
        ctd_005_combined.lon.iloc[9:13],
        ctd_006_combined.lat.iloc[8:12],
        ctd_006_combined.lon.iloc[8:12],
    )
])

grad_t2 = dT2/D2

depth2 = np.linspace(0.15, 1.38,5)


dis2  = np.linspace(0,D2[1],4)*1000

fig, ax1 = plt.subplots(figsize=(10,5))

pcm = ax1.contourf(dis2, depth2, grad_t2,cmap = 'viridis')
contours = ax1.contour(dis2, depth2, grad_t2, colors='black', linewidths=1)
ax1.clabel(contours, inline=True, fontsize=8, fmt="%.2f")  # Add labels to contours
ax1.invert_yaxis()
ax1.set_xlabel("Distance (m)")
ax1.ticklabel_format(style='plain', axis='x', useOffset=False)
ax1.set_ylabel("Depth (m)")
yticks = np.arange(0.15,1.38,0.3)
ax1.set_yticks(yticks)
plt.colorbar(pcm, ax=ax1, label="dT/dy (°C/m)")
plt.title("Horizontal Temperature Gradient Section 2")
plt.show()

#%%


t_31 = ctd_005_combined.loc[17:19,['Temperature (Celsius)']].to_numpy()
t_32 = ctd_006_combined.loc[14:16,['Temperature (Celsius)']].to_numpy()

dT3 = t_31 - t_32

from haversine import haversine
import numpy as np

D3 = np.array([
    haversine((lat1, lon1), (lat2, lon2))
    for lat1, lon1, lat2, lon2 in zip(
        ctd_005_combined.lat.iloc[17:19],
        ctd_005_combined.lon.iloc[17:19],
        ctd_006_combined.lat.iloc[14:16],
        ctd_006_combined.lon.iloc[14:16],
    )
])

grad_t3 = dT3/D3

depth3 = np.linspace(0.15,0.8,3)


dis3  = np.linspace(0,D3[1],2)*1000

fig, ax1 = plt.subplots(figsize=(10,5))

pcm = ax1.contourf(dis3, depth3, grad_t3,cmap = 'viridis')
contours = ax1.contour(dis3, depth3, grad_t3, colors='black', linewidths=1)
ax1.clabel(contours, inline=True, fontsize=8, fmt="%.2f")  # Add labels to contours
ax1.invert_yaxis()
ax1.set_xlabel("Distance (m)")
ax1.ticklabel_format(style='plain', axis='x', useOffset=False)
ax1.set_ylabel("Depth (m)")
yticks = np.linspace(0.15,0.8,6)
ax1.set_yticks(yticks)
plt.colorbar(pcm, ax=ax1, label="dT/dy (°C/m)")
plt.title("Horizontal Temperature Gradient Section 3")
plt.show()

#%%


t_41 = ctd_005_combined.loc[22:28,['Temperature (Celsius)']].to_numpy()
t_42 = ctd_006_combined.loc[17:22,['Temperature (Celsius)']].to_numpy()

min_len = min(len(t_41), len(t_42))
t_41 = t_41[:min_len]
t_42 = t_42[:min_len]

dT4 = t_41 - t_42

from haversine import haversine
import numpy as np

D4 = np.array([
    haversine((lat1, lon1), (lat2, lon2))
    for lat1, lon1, lat2, lon2 in zip(
        ctd_005_combined.lat.iloc[17:21],
        ctd_005_combined.lon.iloc[17:21],
        ctd_006_combined.lat.iloc[17:21],
        ctd_006_combined.lon.iloc[17:21],
    )
])

grad_t4 = dT4/D4

depth4 = np.linspace(0.15, 1.6,6)


dis4  = np.linspace(0,D4[1],4)*1000

fig, ax1 = plt.subplots(figsize=(10,5))

pcm = ax1.contourf(dis4, depth4, grad_t4,cmap = 'viridis')
contours = ax1.contour(dis4, depth4, grad_t4, colors='black', linewidths=1)
ax1.clabel(contours, inline=True, fontsize=8, fmt="%.2f")  # Add labels to contours
ax1.invert_yaxis()
ax1.set_xlabel("Distance (m)")
ax1.ticklabel_format(style='plain', axis='x', useOffset=False)
ax1.set_ylabel("Depth (m)")
yticks = np.linspace(0.15,1.40,6)
ax1.set_yticks(yticks)
plt.colorbar(pcm, ax=ax1, label="dT/dy (°C/m)")
plt.title("Horizontal Temperature Gradient Section 4")
plt.show()
#%%


t_51 = ctd_005_combined.loc[22:28,['Temperature (Celsius)']].to_numpy()
t_52 = ctd_006_combined.loc[17:22,['Temperature (Celsius)']].to_numpy()

min_len = min(len(t_41), len(t_42))
t_51 = t_51[:min_len]
t_52 = t_52[:min_len]

dT5 = t_51 - t_52

from haversine import haversine
import numpy as np

D5 = np.array([
    haversine((lat1, lon1), (lat2, lon2))
    for lat1, lon1, lat2, lon2 in zip(
        ctd_005_combined.lat.iloc[17:21],
        ctd_005_combined.lon.iloc[17:21],
        ctd_006_combined.lat.iloc[17:21],
        ctd_006_combined.lon.iloc[17:21],
    )
])

grad_t5 = dT5/D5

depth5= np.linspace(0.15, 1.6,6)


dis5  = np.linspace(0,D5[1],4)*1000

fig, ax1 = plt.subplots(figsize=(10,5))

pcm = ax1.contourf(dis5, depth5, grad_t5,cmap = 'viridis')
contours = ax1.contour(dis5, depth5, grad_t5, colors='black', linewidths=1)
ax1.clabel(contours, inline=True, fontsize=8, fmt="%.2f")  # Add labels to contours
ax1.invert_yaxis()
ax1.set_xlabel("Distance (m)")
ax1.ticklabel_format(style='plain', axis='x', useOffset=False)
ax1.set_ylabel("Depth (m)")
yticks = np.linspace(0.15,1.40,6)
ax1.set_yticks(yticks)
plt.colorbar(pcm, ax=ax1, label="dT/dy (°C/m)")
plt.title("Horizontal Temperature Gradient Section 4")
plt.show()

#%%
