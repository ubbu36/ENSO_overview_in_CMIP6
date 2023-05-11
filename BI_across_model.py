#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 16:53:09 2021

@author: ullaheede
"""


import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd
import glob as glob
import os
import scipy
import xrscipy.signal as dsp
from pylab import *
import matplotlib.gridspec as gridspec

model_names=['CanESM5','CESM2',\
             'CNRM-CM6-1','HadGEM3-GC31-LL','IPSL-CM6A-LR',\
                    'MIROC6','MIROC-ES2L']

model_names_1=['CanESM5','CESM2',\
             'HadGEM3-GC31-LL',\
                    'MIROC6','MIROC-ES2L']

model_names_2=['CanESM5','CESM2',\
             'CNRM-CM6-1','HadGEM3-GC31-LL','IPSL-CM6A-LR',\
                    'MIROC6']
#model_names=['CanESM5','CESM2-WACCM','CESM2',\
#             'CNRM-CM6-1','CNRM-ESM2-1','FGOALS-g3',\
#                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
#                     'MPI-ESM1-2-LR','UKESM1-0-LL']


nino3_start=180
nino3_stop=280

enso=[0]*len(model_names)
BI=[0]*len(model_names)
ubar=[0]*len(model_names)
wbar=[0]*len(model_names)
vbar=[0]*len(model_names)
alpha_s=[0]*len(model_names)
z_adv=[0]*len(model_names)
ekman=[0]*len(model_names)
thermo=[0]*len(model_names)
noise=[0]*len(model_names_2)

enso_4=[0]*len(model_names_1)
BI_4=[0]*len(model_names_1)
ubar_4=[0]*len(model_names_1)
wbar_4=[0]*len(model_names_1)
alpha_s_4=[0]*len(model_names_1)
z_adv_4=[0]*len(model_names_1)
ekman_4=[0]*len(model_names_1)
thermo_4=[0]*len(model_names_1)
noise_4=[0]*len(model_names_1)

enso_C=[0]*len(model_names)
BI_C=[0]*len(model_names)
ubar_C=[0]*len(model_names)
wbar_C=[0]*len(model_names)
alpha_s_C=[0]*len(model_names)
z_adv_C=[0]*len(model_names)
ekman_C=[0]*len(model_names)
thermo_C=[0]*len(model_names)
noise_C=[0]*len(model_names)


w_nino=180
e_nino=280

L_y=110*1000*10
L_x=110*1000*(e_nino-w_nino)
H=50
rho=4180
cp=1029

#%%

    
for x in range(0,len(model_names)):
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names[x]+'.nc')
    
    control=control['ts']
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    enso_C[x]=low_nino3.std('time')
    
    
  
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1_'+model_names[x]+'.nc')
    control=control['ts'].isel(time=slice(0,85*12))
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    enso1=low_nino3.std('time')

    
    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq2_'+model_names[x]+'.nc')
    control=control['ts'].isel(time=slice(0,85*12))
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    enso2=low_nino3.std('time')


    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq3_'+model_names[x]+'.nc')
    control=control['ts'].isel(time=slice(0,85*12))
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    enso3=low_nino3.std('time')

    enso[x]=(enso1+enso2+enso3)/3
    
    ubar_Cq=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/u_bar_control.nc')
    ubar_C[x]=ubar_Cq['uo']/(L_x)*3600*24*365*(-1)
    ubar1=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/u_bar_ens1.nc')
    ubar2=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/u_bar_ens2.nc')
    ubar3=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/u_bar_ens3.nc')
    ubar[x]=(((ubar1['uo']+ubar2['uo']+ubar3['uo'])/3)/L_x)*3600*24*365*(-1)
 
 #   vbar_Cq=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/v_bar_control.nc')
 #   vbar_C[x]=vbar_Cq['vo']/(L_y)*3600*24*365
 #   vbar1=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/v_bar_ens1.nc')
 #   vbar2=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/v_bar_ens2.nc')
 #   vbar3=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/v_bar_ens3.nc')
 #   vbar[x]=(((vbar1['vo']+vbar2['vo']+vbar3['vo'])/3)/L_y)*3600*24*365

 
    wbar_Cq=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/w_bar_control.nc')
    wbar_C[x]=wbar_Cq['wo']/(H)*3600*24*365*(-1)
    wbar1=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/w_bar_ens1.nc')
    wbar2=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/w_bar_ens2.nc')
    wbar3=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/w_bar_ens3.nc')
    wbar[x]=(((wbar1['wo']+wbar2['wo']+wbar3['wo'])/3)/H)*3600*24*365*(-1)

    alpha_s_C[x]=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/alpha_s_4x_control.npy')/(cp*rho*H)*3600*24*365
    alpha_s1=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/alpha_s_4x_ens1.npy')
    alpha_s2=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/alpha_s_4x_ens2.npy')
    alpha_s3=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/alpha_s_4x_ens3.npy')
    alpha_s[x]=(((alpha_s1+alpha_s2+alpha_s3)/3)/(cp*rho*H))*3600*24*365

    thermoq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/thermo_4x_ens1.npy')
    thermo_C[x]=thermoq[0]*3600*24*365
    thermoq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/thermo_4x_ens1.npy')
    thermo1=thermoq[1]
    thermoq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/thermo_4x_ens2.npy')
    thermo2=thermoq[1]
    thermoq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/thermo_4x_ens3.npy')
    thermo3=thermoq[1]
    thermo[x]=((thermo1+thermo2+thermo3)/3)*3600*24*365
    
    z_advq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/z_adv_4x_ens1.npy')
    z_adv_C[x]=z_advq[0]*3600*24*365
    z_advq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/z_adv_4x_ens1.npy')
    z_adv1=z_advq[1]
    z_advq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/z_adv_4x_ens2.npy')
    z_adv2=z_advq[1]
    z_advq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/z_adv_4x_ens3.npy')
    z_adv3=z_advq[1]
    z_adv[x]=((z_adv1+z_adv2+z_adv3)/3)*3600*24*365
    
    ekmanq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/ekman_4x_ens1.npy')
    ekman_C[x]=ekmanq[0]*3600*24*365
    ekmanq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/ekman_4x_ens1.npy')
    ekman1=ekmanq[1]
    ekmanq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/ekman_4x_ens2.npy')
    ekman2=ekmanq[1]
    ekmanq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/ekman_4x_ens3.npy')
    ekman3=ekmanq[1]
    ekman[x]=((ekman1+ekman2+ekman3)/3)*3600*24*365
    
    BI_C[x]=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/BI_control.npy')
    BI1=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/BI_4x_ens1.npy')
    BI2=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/BI_4x_ens2.npy')
    BI3=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names[x]+'/BI_4x_ens3.npy')
    BI[x]=((BI1+BI2+BI3)/3)
    

    noise_t=xr.open_dataset('/Users/ullaheede_1/Downloads/calc_noise_'+model_names[x]+'_control.nc')
    noise_C[x]=noise_t['__xarray_dataarray_variable__']
#%%
for x in range(0,len(model_names_2)):
    noise_1=xr.open_dataset('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_2[x]+'_1.nc')
    noise_1=noise_1['__xarray_dataarray_variable__']
    noise_2=xr.open_dataset('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_2[x]+'_2.nc')
    noise_2=noise_2['__xarray_dataarray_variable__']
    noise_3=xr.open_dataset('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_2[x]+'_3.nc')
    noise_3=noise_3['__xarray_dataarray_variable__']
    noise[x]=(noise_1+noise_2+noise_3)/3
#    enso_C[x]=ensoC
    
for x in range(0,len(model_names_1)):
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_1[x]+'.nc')
    control=control['ts']
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    ensoC=low_nino3.std('time')
#    enso_C[x]=ensoC
    
  
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq4_'+model_names_1[x]+'.nc')
    control=control['ts'].isel(time=slice(0,150*12))
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    enso1=low_nino3.std('time')



    enso_4[x]=(enso1)
    
    ubar_Cq=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/u_bar_control.nc')
#    ubar_C[x]=ubar_Cq['uo']
    ubar1=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/u_bar_4x.nc')

    ubar_4[x]=((ubar1['uo'])/L_x)*3600*24*365*(-1)
    
    wbar_Cq=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/w_bar_control.nc')
#    wbar_C[x]=wbar_Cq['wo']
    wbar1=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/w_bar_4x.nc')

    wbar_4[x]=((wbar1['wo'])/H)*3600*24*365*(-1)

 #  alpha_s_C[x]=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/alpha_s_4x_control.npy')
    alpha_s1=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/alpha_s_4x_4x.npy')
    alpha_s_4[x]=((alpha_s1)/(cp*rho*H))*3600*24*365

    thermoq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/thermo_4x_ens1.npy')
#    thermo_C[x]=thermoq[0]
    thermoq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/thermo_4x_4x.npy')
    thermo1=thermoq[1]

    thermo_4[x]=(thermo1)*3600*24*365
    
    z_advq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/z_adv_4x_ens1.npy')
#    z_adv_C[x]=z_advq[0]
    z_advq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/z_adv_4x_4x.npy')
    z_adv1=z_advq[1]

    z_adv_4[x]=(z_adv1)*3600*24*365
    
    ekmanq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/ekman_4x_ens1.npy')
#    ekman_C[x]=ekmanq[0]
    ekmanq=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/ekman_4x_4x.npy')
    ekman1=ekmanq[1]

    ekman_4[x]=(ekman1)*3600*24*365
    
 #   BI_C[x]=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/BI_control.npy')
    BI1=np.load('/Users/ullaheede_1/Documents/ENSO_project/'+model_names_1[x]+'/BI_4x_4x.npy')

    BI_4[x]=(BI1)
    

    noise_t=xr.open_dataset('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_1[x]+'_4x.nc')
    noise_4[x]=noise_t['__xarray_dataarray_variable__']

#%%
plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'hatch.color': '0.1'})  
x = np.arange(len(model_names))  # the label locations
x1=array([ 0, 1,3,5,6])
x2=array([ 0,1,2,3,4,5])

width = 0.2  # the width of the bars
fig = figure(figsize=(47,38))
gs = gridspec.GridSpec(4, 5)
ax1 = plt.subplot(gs[0, 0:1])
ax2 = plt.subplot(gs[0, 1:2])
ax3 = plt.subplot(gs[0, 2:3])
ax4 = plt.subplot(gs[0, 3:4])
ax5 = plt.subplot(gs[1, 0:1])
ax6 = plt.subplot(gs[1, 1:2])
ax7 = plt.subplot(gs[1, 2:3])
ax8 = plt.subplot(gs[1, 3:4])

ax9 = plt.subplot(gs[0, 4:5])

fig = gcf()
gs.tight_layout(fig,h_pad=10,w_pad=1)
ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7,ax8,ax9]

plt.figtext(0.032, 0.982, 'a)')
plt.figtext(0.228, 0.982, 'b)')
plt.figtext(0.43, 0.982, 'c)')
plt.figtext(0.63, 0.982, 'd)')
plt.figtext(0.82, 0.982, 'e)')

plt.figtext(0.032, 0.7, 'f)')
plt.figtext(0.228, 0.7, 'g)')
plt.figtext(0.43, 0.7, 'h)')
plt.figtext(0.63, 0.7, 'i)')

ax[0].bar(x - width, enso_C, width,label='PiControl',color='k')
ax[0].bar(x1, enso_4, width,label='abrupt-4xCO2',color='#1f77b4')
ax[0].bar(x+width, enso, width,label='ssp585',color='#d62728')

ax[0].set_xticks(range(len(model_names)))

ax[0].set_xticklabels(model_names,rotation='vertical')

ax[0].set_title('ENSO SST',fontsize=45)

#ax[0].set_ylabel('eq. std')
#ax[0].legend(ncol=2,fontsize=38)

ax[1].bar(x - width, ubar_C, width,label='PiControl',color='k')
ax[1].bar(x+width, ubar, width,label='ssp585',color='#d62728')
ax[1].bar(x1, ubar_4, width,label='abrupt-4xCO2',color='#1f77b4')

ax[1].set_xticks(range(len(model_names)))

ax[1].set_xticklabels(model_names,rotation='vertical')

ax[1].set_title('$\overline{u}$/L$_x$',fontsize=45)

#ax[1].set_ylabel('eq. std')
#ax[0].legend(ncol=2,fontsize=38)

ax[2].bar(x - width, wbar_C, width,label='PiControl',color='k')
ax[2].bar(x+width, wbar, width,label='ssp585',color='#d62728')
ax[2].bar(x1, wbar_4, width,label='abrupt-4xCO2',color='#1f77b4')

ax[2].set_xticks(range(len(model_names)))

ax[2].set_xticklabels(model_names,rotation='vertical')

ax[2].set_title('$\overline{w}$/H$_m$',fontsize=45)

#ax[2].set_ylabel('eq. std')
#ax[0].legend(ncol=2,fontsize=38)

ax[3].bar(x - width, alpha_s_C, width,label='PiControl',color='k')
ax[3].bar(x+width, alpha_s, width,label='ssp585',color='#d62728')
ax[3].bar(x1, alpha_s_4, width,label='abrupt-4xCO2',color='#1f77b4')

ax[3].set_xticks(range(len(model_names)))

ax[3].set_xticklabels(model_names,rotation='vertical')

ax[3].set_title('$\\alpha$',fontsize=45)

#ax[3].set_ylabel('eq. std')
#ax[0].legend(ncol=2,fontsize=38)



ax[4].bar(x - width, BI_C, width,label='PiControl',color='k')
ax[4].bar(x+width, BI, width,label='ssp585',color='#d62728')
ax[4].bar(x1, BI_4, width,label='abrupt-4xCO2',color='#1f77b4')

ax[4].set_xticks(range(len(model_names)))

ax[4].set_xticklabels(model_names,rotation='vertical')

ax[4].set_title('Bjerknes Index',fontsize=45)

#ax[4].set_ylabel('eq. std')


ax[5].bar(x - width, z_adv_C, width,label='PiControl',color='k')
ax[5].bar(x+width, z_adv, width,label='ssp585',color='#d62728')
ax[5].bar(x1, z_adv_4, width,label='abrupt-4xCO2',color='#1f77b4')

ax[5].set_xticks(range(len(model_names)))

ax[5].set_xticklabels(model_names,rotation='vertical')

ax[5].set_title('zonal adv. feedback',fontsize=45)

#ax[5].set_ylabel('eq. std')

ax[6].bar(x - width, ekman_C, width,label='PiControl',color='k')
ax[6].bar(x+width, ekman, width,label='ssp585',color='#d62728')
ax[6].bar(x1, ekman_4, width,label='abrupt-4xCO2',color='#1f77b4')

ax[6].set_xticks(range(len(model_names)))

ax[6].set_xticklabels(model_names,rotation='vertical')

ax[6].set_title('ekman feedback',fontsize=45)

#ax[6].set_ylabel('eq. std')

ax[7].bar(x - width, thermo_C, width,label='PiControl',color='k')
ax[7].bar(x+width, thermo, width,label='ssp585',color='#d62728')
ax[7].bar(x1, thermo_4, width,label='abrupt-4xCO2',color='#1f77b4')

ax[7].set_xticks(range(len(model_names)))

ax[7].set_xticklabels(model_names,rotation='vertical')

ax[7].set_title('thermocl. feedback',fontsize=45)

#ax[7].set_ylabel('eq. std')


#ax[8].set_ylabel('eq. std')

ax[8].bar(x - width, noise_C, width,label='PiControl',color='k')
ax[8].bar(x2+width, noise, width,label='ssp585',color='#d62728')
ax[8].bar(x1, noise_4, width,label='abrupt-4xCO2',color='#1f77b4')

ax[8].set_xticks(range(len(model_names)))

ax[8].set_xticklabels(model_names,rotation='vertical')

ax[8].set_title('noise',fontsize=45)

#ax[8].set_ylabel('eq. std')
#ax[0].legend(ncol=2,fontsize=38)

#ax[4].legend(ncol=2,fontsize=22)

lines = []
labels = []


axLine, axLabel = ax1.get_legend_handles_labels()
lines.extend(axLine)
labels.extend(axLabel)

    
fig.legend(lines, labels,           
           loc = (0.85,0.36))

plt.show()
#%%
alpha_s_C=np.array(alpha_s_C)
alpha_s=np.array(alpha_s)

enso_C=np.array(enso_C)
enso=np.array(enso)
