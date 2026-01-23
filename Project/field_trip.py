#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 10:56:58 2026

@author: bobco-08
"""

# import numpy as np
# import xarray as xr
# import matplotlib.pyplot as plt
# import pandas as pd

# file_1 = '/home/bobco-08/24cl05012/field trip/ctd_profiles_combined.csv'
# file_2 = '/home/bobco-08/24cl05012/field trip/ctd_profiles_combined-1.csv'

# data_1 =  pd.read_csv(file_1)

# data_2 = pd.read_csv(file_2)

# #%%

# para_005 = data_1.loc[0:196,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)','lat','lon','cast_time_local']]


# para_1_005 = para_005.loc[0:8,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
# para_12_005 = para_005.loc[9:17,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# loc_1 = para_005.loc[0:8,['cast_time_local','lat','lon']]

# avg_para_1_005 = (para_1_005.reset_index(drop=True)+ para_12_005.reset_index(drop=True))/2

# avg_1_005 = pd.concat([loc_1,avg_para_1_005,],axis =1)

# #%%
# para_2_005 = para_005.loc[18:25,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
# para_22_005 = para_005.loc[26:33,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# loc_2 = para_005.loc[18:25,['cast_time_local','lat','lon']].reset_index(drop =True)

   
# para_23_005 = para_005.loc[26:33,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_24_005 = para_005.loc[34:41,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_25_005 = para_005.loc[34:41,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_26_005 = para_005.loc[50:57,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_27_005 = para_005.loc[58:65,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]



# avg_para_2_005 = (para_2_005.reset_index(drop=True)+ para_22_005.reset_index(drop=True)+para_23_005.reset_index(drop=True)+para_24_005.reset_index(drop=True)+para_25_005.reset_index(drop=True)+para_26_005.reset_index(drop=True)+para_27_005.reset_index(drop=True))/7

# avg_2_005 = pd.concat([loc_2,avg_para_2_005,],axis =1)

# #%%

# para_3_005 = para_005.loc[66:70,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
# para_32_005 = para_005.loc[71:75,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
# para_33_005 = para_005.loc[76:80,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# loc_3 = para_005.loc[66:70,['cast_time_local','lat','lon']].reset_index(drop =True)


# avg_para_3_005 = (para_3_005.reset_index(drop=True)+ para_32_005.reset_index(drop=True)+para_33_005.reset_index(drop=True))/3

# avg_3_005 = pd.concat([loc_3,avg_para_3_005,],axis =1)


# #%%

# para_4_005 = para_005.loc[81:87,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
# para_42_005 = para_005.loc[88:94,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

   
# para_43_005 = para_005.loc[95:101,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_44_005 = para_005.loc[102:108,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_45_005 = para_005.loc[109:115,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]


# loc_4 = para_005.loc[81:87,['cast_time_local','lat','lon']].reset_index(drop =True)


# avg_para_4_005 = (para_4_005.reset_index(drop=True)+ para_42_005.reset_index(drop=True)+para_43_005.reset_index(drop=True)+para_44_005.reset_index(drop=True)+para_45_005.reset_index(drop=True))/5

# avg_4_005 = pd.concat([loc_4,avg_para_4_005,],axis =1)

# #%%

# para_5_005 = para_005.loc[116:123,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_52_005 = para_005.loc[124:131,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_53_005 = para_005.loc[132:139,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]


# loc_5 = para_005.loc[116:123,['cast_time_local','lat','lon']].reset_index(drop =True)


# avg_para_5_005 = (para_5_005.reset_index(drop=True)+ para_52_005.reset_index(drop=True)+para_53_005.reset_index(drop=True))/3

# avg_5_005 = pd.concat([loc_5,avg_para_5_005,],axis =1)


# #%%

# para_6_005 = para_005.loc[140:146,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_62_005 = para_005.loc[147:153,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_63_005 = para_005.loc[154:160,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# loc_6 = para_005.loc[140:146,['cast_time_local','lat','lon']].reset_index(drop =True)

# avg_para_6_005 = (para_6_005.reset_index(drop=True)+ para_62_005.reset_index(drop=True)+para_63_005.reset_index(drop=True))/3

# avg_6_005 = pd.concat([loc_6,avg_para_6_005,],axis =1)

# #%%

# para_7_005 = para_005.loc[161:167,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_72_005 = para_005.loc[168:174,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_73_005 = para_005.loc[175:181,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# loc_7 = para_005.loc[161:167,['cast_time_local','lat','lon']].reset_index(drop =True)

# avg_para_7_005 = (para_7_005.reset_index(drop=True)+ para_72_005.reset_index(drop=True)+para_73_005.reset_index(drop=True))/3

# avg_7_005 = pd.concat([loc_7,avg_para_7_005,],axis =1)

# #%%

# para_8_005 = para_005.loc[182:186,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_82_005 = para_005.loc[187:191,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# para_83_005 = para_005.loc[192:196,['Depth (Meter)','Temperature (Celsius)','Salinity (Practical Salinity Scale)','Sound velocity (Meters per Second)','Density (Kilograms per Cubic Meter)']]

# loc_8 = para_005.loc[182:186,['cast_time_local','lat','lon']].reset_index(drop =True)

# avg_para_8_005 = (para_8_005.reset_index(drop=True)+ para_82_005.reset_index(drop=True)+para_83_005.reset_index(drop=True))/3

# avg_8_005 = pd.concat([loc_8,avg_para_8_005,],axis =1)

# #%%

# ctd_005_combined =  pd.concat([avg_1_005,avg_2_005,avg_3_005,avg_4_005,avg_5_005,avg_6_005,avg_7_005,avg_8_005]).reset_index(drop=True)

# ctd_005_combined.to_csv('/home/bobco-08/24cl05012/field trip/ctd_005_results.csv', index=False)

#%% temperature profile
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

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
axs[0,0].set_title('19.664N , 85.2138E (2025-12-09 14:50:27)', pad=18)
axs[0,0].set_xlabel('Temp (°C)')
axs[0,0].set_ylabel('Depth (m)')
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(region_2['Temperature (Celsius)'], region_2['Depth (Meter)'], 'ro-')
axs[0,1].set_title('19.6903N , 85.2101E (2025-12-09 16:14:38)', pad=18)
axs[0,1].set_xlabel('Temp (°C)')
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(region_3['Temperature (Celsius)'], region_3['Depth (Meter)'], 'ro-')
axs[0,2].set_title('19.7159N , 85.1956E (2025-12-09 16:56:03)', pad=18)
axs[0,2].set_xlabel('Temp (°C)')
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(region_4['Temperature (Celsius)'], region_4['Depth (Meter)'], 'ro-')
axs[0,3].set_title('19.7266N , 85.2435E (2025-12-10 09:33:05)', pad=18)
axs[0,3].set_xlabel('Temp (°C)')
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(region_5['Temperature (Celsius)'], region_5['Depth (Meter)'], 'ro-')
axs[1,0].set_title('19.729N , 85.2755E (2025-12-10 10:18:35)', pad=18)
axs[1,0].set_xlabel('Temp (°C)')
axs[1,0].set_ylabel('Depth (m)')
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(region_6['Temperature (Celsius)'], region_6['Depth (Meter)'], 'ro-')
axs[1,1].set_title('19.7276N , 85.3145E (2025-12-10 11:18:26)', pad=18)
axs[1,1].set_xlabel('Temp (°C)')
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(region_7['Temperature (Celsius)'], region_7['Depth (Meter)'], 'ro-')
axs[1,2].set_title('19.7234N , 85.2187E (2025-12-10 12:48:51)', pad=18)
axs[1,2].set_xlabel('Temp (°C)')
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(region_8['Temperature (Celsius)'], region_8['Depth (Meter)'], 'ro-')
axs[1,3].set_title('19.7161N , 85.196E (2025-12-10 13:41:05)', pad=18)
axs[1,3].set_xlabel('Temp (°C)')
axs[1,3].grid(True, alpha=0.3)

# ---------- Common axis formatting ----------
temp_ticks = np.linspace(21.8, 23.8, 9)
dep_ticks = np.linspace(0.15,2.8,5)
for ax in axs.flat:
    ax.set_xlim(21.8, 23.8)
    ax.set_ylim(2.8, 0.15)
    ax.set_xticks(temp_ticks)
    ax.set_yticks(dep_ticks )
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
plt.tight_layout()
plt.show()



#%% salinity

fig, axs = plt.subplots(2, 4, figsize=(20,15), sharey=True)

# -------- Region 1 --------
axs[0,0].plot(region_1['Salinity (Practical Salinity Scale)'], region_1['Depth (Meter)'], 'bo-',)
axs[0,0].set_title('19.664N , 85.2138E (2025-12-09 14:50:27)', pad=18)
axs[0,0].set_xlabel('PSU')
axs[0,0].set_ylabel('Depth (m)')
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(region_2['Salinity (Practical Salinity Scale)'], region_2['Depth (Meter)'], 'bo-')
axs[0,1].set_title('19.6903N , 85.2101E (2025-12-09 16:14:38)', pad=18)
axs[0,1].set_xlabel('PSU')
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(region_3['Salinity (Practical Salinity Scale)'], region_3['Depth (Meter)'], 'bo-')
axs[0,2].set_title('19.7159N , 85.1956E (2025-12-09 16:56:03)', pad=18)
axs[0,2].set_xlabel('PSU')
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(region_4['Salinity (Practical Salinity Scale)'], region_4['Depth (Meter)'], 'bo-')
axs[0,3].set_title('19.7266N , 85.2435E (2025-12-10 09:33:05)', pad=18)
axs[0,3].set_xlabel('PSU')
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(region_5['Salinity (Practical Salinity Scale)'], region_5['Depth (Meter)'], 'bo-')
axs[1,0].set_title('19.729N , 85.2755E (2025-12-10 10:18:35)', pad=18)
axs[1,0].set_xlabel('PSU')
axs[1,0].set_ylabel('Depth (m)')
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(region_6['Salinity (Practical Salinity Scale)'], region_6['Depth (Meter)'], 'bo-')
axs[1,1].set_title('19.7276N , 85.3145E (2025-12-10 11:18:26)', pad=18)
axs[1,1].set_xlabel('Temp (°C)')
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(region_7['Salinity (Practical Salinity Scale)'], region_7['Depth (Meter)'], 'bo-')
axs[1,2].set_title('19.7234N , 85.2187E (2025-12-10 12:48:51)', pad=18)
axs[1,2].set_xlabel('PSU')
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(region_8['Salinity (Practical Salinity Scale)'], region_8['Depth (Meter)'], 'bo-')
axs[1,3].set_title('19.7161N , 85.196E (2025-12-10 13:41:05)', pad=18)
axs[1,3].set_xlabel('PSU')
axs[1,3].grid(True, alpha=0.3)

#---------- Common axis formatting ----------
sal_ticks = np.linspace(7.8, 9.5, 9)
dep_ticks = np.linspace(0.15,2.8,5)
for ax in axs.flat:
    ax.set_xlim(7.8, 9.5)
    ax.set_ylim(2.8, 0.15)
    ax.set_xticks(sal_ticks)
    ax.set_yticks(dep_ticks )
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
plt.tight_layout()
plt.show()

#%%

import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 4, figsize=(20,15), sharey=True)

# -------- Region 1 --------
axs[0,0].plot(region_1['Sound velocity (Meters per Second)'], region_1['Depth (Meter)'],'o-')
axs[0,0].set_title('19.664N , 85.2138E (2025-12-09 14:50:27)', pad=18)
axs[0,0].set_xlabel('Sound velocity (m/s)')
axs[0,0].set_ylabel('Depth (m)')
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(region_2['Sound velocity (Meters per Second)'], region_2['Depth (Meter)'],'o-')
axs[0,1].set_title('19.6903N , 85.2101E (2025-12-09 16:14:38)', pad=18)
axs[0,1].set_xlabel('Sound velocity (m/s)')
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(region_3['Sound velocity (Meters per Second)'], region_3['Depth (Meter)'], 'o-')
axs[0,2].set_title('19.7159N , 85.1956E (2025-12-09 16:56:03)', pad=18)
axs[0,2].set_xlabel('Sound velocity (m/s)')
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(region_4['Sound velocity (Meters per Second)'], region_4['Depth (Meter)'], 'o-')
axs[0,3].set_title('19.7266N , 85.2435E (2025-12-10 09:33:05)', pad=18)
axs[0,3].set_xlabel('Sound velocity (m/s)')
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(region_5['Sound velocity (Meters per Second)'], region_5['Depth (Meter)'], 'o-')
axs[1,0].set_title('19.729N , 85.2755E (2025-12-10 10:18:35)', pad=18)
axs[1,0].set_xlabel('Sound velocity (m/s)')
axs[1,0].set_ylabel('Depth (m)')
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(region_6['Sound velocity (Meters per Second)'], region_6['Depth (Meter)'], 'o-')
axs[1,1].set_title('19.7276N , 85.3145E (2025-12-10 11:18:26)', pad=18)
axs[1,1].set_xlabel('Sound velocity (m/s)')
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(region_7['Sound velocity (Meters per Second)'], region_7['Depth (Meter)'], 'o-')
axs[1,2].set_title('19.7234N , 85.2187E (2025-12-10 12:48:51)', pad=18)
axs[1,2].set_xlabel('Sound velocity (m/s)')
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(region_8['Sound velocity (Meters per Second)'], region_8['Depth (Meter)'], 'o-')
axs[1,3].set_title('19.7161N , 85.196E (2025-12-10 13:41:05)', pad=18)
axs[1,3].set_xlabel('Sound velocity (m/s)')
axs[1,3].grid(True, alpha=0.3)

#---------- Common axis formatting ----------
# You can adjust x-limits according to your sound velocity range
sound_ticks = np.linspace(1496, 1504, 9)  # example: typical ocean sound speed range
dep_ticks = np.linspace(0.15,2.8,5)
for ax in axs.flat:
    ax.set_xlim(1496, 1504)
    ax.set_ylim(2.8, 0.15)  # invert depth axis
    ax.set_xticks(sound_ticks)
    ax.set_yticks(dep_ticks)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

plt.tight_layout()
plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 4, figsize=(20,15), sharey=True)

# -------- Region 1 --------
axs[0,0].plot(region_1['Density (Kilograms per Cubic Meter)'], region_1['Depth (Meter)'], 'ko-')
axs[0,0].set_title('19.664N , 85.2138E (2025-12-09 14:50:27)', pad=18)
axs[0,0].set_xlabel('Density (kg/m³)')
axs[0,0].set_ylabel('Depth (m)')
axs[0,0].grid(True, alpha=0.3)

# -------- Region 2 --------
axs[0,1].plot(region_2['Density (Kilograms per Cubic Meter)'], region_2['Depth (Meter)'], 'ko-')
axs[0,1].set_title('19.6903N , 85.2101E (2025-12-09 16:14:38)', pad=18)
axs[0,1].set_xlabel('Density (kg/m³)')
axs[0,1].grid(True, alpha=0.3)

# -------- Region 3 --------
axs[0,2].plot(region_3['Density (Kilograms per Cubic Meter)'], region_3['Depth (Meter)'], 'ko-')
axs[0,2].set_title('19.7159N , 85.1956E (2025-12-09 16:56:03)', pad=18)
axs[0,2].set_xlabel('Density (kg/m³)')
axs[0,2].ticklabel_format(style='plain', axis='x', useOffset=False)
axs[0,2].grid(True, alpha=0.3)

# -------- Region 4 --------
axs[0,3].plot(region_4['Density (Kilograms per Cubic Meter)'], region_4['Depth (Meter)'], 'ko-')
axs[0,3].set_title('19.7266N , 85.2435E (2025-12-10 09:33:05)', pad=18)
axs[0,3].set_xlabel('Density (kg/m³)')
axs[0,3].grid(True, alpha=0.3)

# -------- Region 5 --------
axs[1,0].plot(region_5['Density (Kilograms per Cubic Meter)'], region_5['Depth (Meter)'], 'ko-')
axs[1,0].set_title('19.729N , 85.2755E (2025-12-10 10:18:35)', pad=18)
axs[1,0].set_xlabel('Density (kg/m³)')
axs[1,0].set_ylabel('Depth (m)')
axs[1,0].grid(True, alpha=0.3)

# -------- Region 6 --------
axs[1,1].plot(region_6['Density (Kilograms per Cubic Meter)'], region_6['Depth (Meter)'], 'ko-')
axs[1,1].set_title('19.7276N , 85.3145E (2025-12-10 11:18:26)', pad=18)
axs[1,1].set_xlabel('Density (kg/m³)')
axs[1,1].grid(True, alpha=0.3)

# -------- Region 7 --------
axs[1,2].plot(region_7['Density (Kilograms per Cubic Meter)'], region_7['Depth (Meter)'], 'ko-')
axs[1,2].set_title('19.7234N , 85.2187E (2025-12-10 12:48:51)', pad=18)
axs[1,2].set_xlabel('Density (kg/m³)')
axs[1,2].grid(True, alpha=0.3)

# -------- Region 8 --------
axs[1,3].plot(region_8['Density (Kilograms per Cubic Meter)'], region_8['Depth (Meter)'], 'ko-')
axs[1,3].set_title('19.7161N , 85.196E (2025-12-10 13:41:05)', pad=18)
axs[1,3].set_xlabel('Density (kg/m³)')
axs[1,3].grid(True, alpha=0.3)

#---------- Common axis formatting ----------
# You can adjust x-limits according to your density range
density_ticks = np.arange(1003, 1005.5,0.5)  # example typical seawater density
dep_ticks = np.linspace(0.15,2.8,5)

for ax in axs.flat:
    ax.set_xlim(1003, 1005)
    ax.set_ylim(2.8, 0.15)  # invert depth axis
    ax.set_xticks(density_ticks)
    ax.set_yticks(dep_ticks)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

plt.tight_layout()
plt.show()

#%%


import pandas as pd

# List of all regions
regions = [region_1, region_2, region_3, region_4,
           region_5, region_6, region_7, region_8]

# Optional: names for each region
region_names = [f'Region {i}' for i in range(1, 9)]

# Initialize variables
max_sal = -float('inf')
min_sal = float('inf')
max_region = ''
min_region = ''

# Loop through each region
for reg, name in zip(regions, region_names):
    sal_max = reg['Density (Kilograms per Cubic Meter)'].max()
    sal_min = reg['Density (Kilograms per Cubic Meter)'].min()
    
    if sal_max > max_sal:
        max_sal = sal_max
        max_region = name
        
    if sal_min < min_sal:
        min_sal = sal_min
        min_region = name

print(f'Maximum salinity: {max_sal} in {max_region}')
print(f'Minimum salinity: {min_sal} in {min_region}')

#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

ctd_005_combined = pd.read_csv('/home/bobco-08/24cl05012/field trip/ctd_005_results.csv')

# -------- Regions --------
regions = [
    ctd_005_combined.iloc[0:9],
    ctd_005_combined.iloc[9:17],
    ctd_005_combined.iloc[17:22],
    ctd_005_combined.iloc[22:29],
    ctd_005_combined.iloc[29:37],
    ctd_005_combined.iloc[37:44],
    ctd_005_combined.iloc[44:51],
    ctd_005_combined.iloc[51:56],
]

labels = [
    '19.664N , 85.2138E (2025-12-09 14:50:27)',
    '19.6903N , 85.2101E (2025-12-09 16:14:38)',
    '19.7159N , 85.1956E (2025-12-09 16:56:03)',
    '19.7266N , 85.2435E (2025-12-10 09:33:05)',
    '19.729N , 85.2755E (2025-12-10 10:18:35)',
    '19.7276N , 85.3145E (2025-12-10 11:18:26)',
    '19.7234N , 85.2187E (2025-12-10 12:48:51)',
    '19.7161N , 85.196E (2025-12-10 13:41:05)'
]

# -------- Single Plot --------
fig, ax = plt.subplots(figsize=(8, 10))

for reg, lab in zip(regions, labels):
    ax.plot(reg['Temperature (Celsius)'],
            reg['Depth (Meter)'],
            marker='o',
            label=lab)

# -------- Common Axis Formatting --------
temp_ticks = np.linspace(21.8, 23.8, 9)
dep_ticks = np.linspace(0.15, 3, 9)

ax.set_xlim(21.8, 23.8)
ax.set_ylim(3,0.15)          # invert depth
ax.set_xticks(temp_ticks)
ax.set_yticks(dep_ticks)

ax.set_xlabel('Temperature (°C)')
ax.set_ylabel('Depth (m)')
ax.set_title('CTD Temperature Profiles – All Stations')
ax.grid(True, alpha=0.3)

ax.legend(fontsize=9)
plt.tight_layout()
plt.show()

#%%

import matplotlib.pyplot as plt
import numpy as np

# Common ticks
temp_ticks = np.linspace(21.8, 23.8, 9)
sal_ticks = np.linspace(7.8, 9.5, 9)
dens_ticks = np.arange(1003, 1005.5,0.5)
dep_ticks = np.linspace(0.15, 2.8, 5)

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

fig, axs = plt.subplots(2, 4, figsize=(24, 12),dpi = 200, sharey=True)

for i, (region, ax) in enumerate(zip(regions, axs.flat)):
    # 1. PRIMARY X-AXIS: Salinity (Bottom)
    ax.plot(region['Salinity (Practical Salinity Scale)'], region['Depth (Meter)'],
            'bo-', label='Salinity')
    ax.set_xlabel('Salinity (PSU)', color='blue')
    ax.set_xlim(sal_ticks[0], sal_ticks[-1])
    ax.set_ylim(dep_ticks[-1], dep_ticks[0])  # Invert depth
    ax.tick_params(axis='x', labelcolor='blue')

    # 2. SECOND X-AXIS: Temperature (Top)
    ax_temp = ax.twiny()
    ax_temp.plot(region['Temperature (Celsius)'], region['Depth (Meter)'],
                 'ro-', label='Temperature')
    ax_temp.set_xlabel('Temperature (°C)', color='red')
    ax_temp.set_xlim(temp_ticks[0], temp_ticks[-1])
    ax_temp.tick_params(axis='x', labelcolor='red')

    # 3. THIRD X-AXIS: Density (Offset Top or Bottom)
    ax_dens = ax.twiny()
    # Offset the third axis position so it's above the Temperature axis
    ax_dens.spines['top'].set_position(('outward', 40)) 
    ax_dens.plot(region['Density (Kilograms per Cubic Meter)'], region['Depth (Meter)'],
                 'go-', label='Density')
    ax_dens.set_xlabel('Density (kg/m³)', color='green')
    ax_dens.set_xlim(dens_ticks[0], dens_ticks[-1])
    ax_dens.tick_params(axis='x', labelcolor='green')

    # General Formatting
    ax.set_ylabel('Depth (Meter)')
    ax.grid(True, alpha=0.3)
    ax.set_title(titles[i], pad=10) # Increased pad to accommodate offset axis

plt.subplots_adjust(
    left=0.05, 
    right=0.95, 
    top=0.85, 
    bottom=0.1, 
    hspace=0.7,  # Increase this for vertical gap between rows
    wspace=0.3   # Increase this for horizontal gap between columns
)
plt.show()