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

model_names=['ACCESS-ESM1-5','CESM2',\
             'IPSL-CM6A-LR']


    
#model_names=['CanESM5','CESM2-WACCM','CESM2',\
#             'CNRM-CM6-1','CNRM-ESM2-1','FGOALS-g3',\
#                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
#                     'MPI-ESM1-2-LR','UKESM1-0-LL']


nino3_start=200
nino3_stop=280

nino4_start=160
nino4_stop=210
    
e1=200
e2=280

w1=120
w2=180
x=20


enso_control=[0]*len(model_names)
enso_ssp585_0=[0]*len(model_names)
enso_ssp585_1=[0]*len(model_names)
enso_ssp585_2=[0]*len(model_names)
enso_ssp585_3=[0]*len(model_names)
enso_ssp585_4=[0]*len(model_names)

noise_control=[0]*len(model_names)
noise_ssp585_0=[0]*len(model_names)
noise_ssp585_1=[0]*len(model_names)
noise_ssp585_2=[0]*len(model_names)
noise_ssp585_3=[0]*len(model_names)
noise_ssp585_4=[0]*len(model_names)

#%%

    
for x in range(0,len(model_names)):
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names[x]+'.nc')
    control=control['ts']
    control=control.assign_coords(time=list(range(len(control.time))))
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=control.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    ensoC=low_nino3.std('time')
    
    
    enso_control[x]=ensoC
    
  
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq_4x_200yr_'+model_names[x]+'.nc')
    control=control['ts'].isel(time=slice(0,200*12))
    control=control.assign_coords(time=list(range(len(control.time))))
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=control.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')
    
    
    enso_ssp585_0[x]=enso1
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq_4x_400yr_'+model_names[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,200*12))
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_1[x]=enso1
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq_4x_600yr_'+model_names[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,200*12))
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso2=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_2[x]=enso2
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq_4x_800yr_'+model_names[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,200*12))
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    
    enso3=low_nino3.std('time')
    
    enso_ssp585_3[x]=enso3
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq_4x_1000yr_'+model_names[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,200*12))
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso4=low_nino3.std('time')

    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_4[x]=enso4
    test=xr.concat([enso1,enso2,enso3], dim='ens')



    plt.plot(ensoC)
    plt.plot(enso_ssp)
#plt.plot(control.mean('lat').mean('time'))
#plt.plot(ens1.mean('lat').mean('time'))

#%%


    
for x in range(0,len(model_names)-1):
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names[x]+'_control.nc')
    control=control['__xarray_dataarray_variable__'].std('time')
    control=control.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')

    ensoC=control
    
    
    noise_control[x]=ensoC
    
  
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names[x]+'_4x_200yr.nc')
    control=control['__xarray_dataarray_variable__'].std('time')
    control=control.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')

    ensoC=control
    
    
    noise_ssp585_0[x]=ensoC
    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names[x]+'_4x_400yr.nc')
    control=control['__xarray_dataarray_variable__'].std('time')
    control=control.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')

    ensoC=control
    
    
    noise_ssp585_1[x]=ensoC
    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names[x]+'_4x_600yr.nc')

    control=control['__xarray_dataarray_variable__'].std('time')
    control=control.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')

    ensoC=control
    
    
    noise_ssp585_2[x]=ensoC
    

    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names[x]+'_4x_800yr.nc')

    control=control['__xarray_dataarray_variable__'].std('time')
    control=control.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')

    ensoC=control
    
    
    noise_ssp585_3[x]=ensoC
    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names[x]+'_4x_1000yr.nc')
    control=control['__xarray_dataarray_variable__'].std('time')
    control=control.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')

    ensoC=control
    
    
    noise_ssp585_4[x]=ensoC
 #   test=xr.concat([enso1,enso2,enso3], dim='ens')



    plt.plot(ensoC)
    plt.plot(enso_ssp)
#%%
plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'hatch.color': '0.1'})  
x = np.arange(len(model_names))
x1 = np.arange(len(model_names)-1)   # the label locations
#x1=array([ 0,  1, 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,16,17,18,19])
#x2=array([ 0,  1, 2, 3,  5,  6,  7,  8,  9, 10,11, 14, 15,16,17,18,19])
width = 0.15  # the width of the bars
fig = figure(figsize=(47,25))
gs = gridspec.GridSpec(2, 3)
ax1 = plt.subplot(gs[0, 0:3])
ax2 = plt.subplot(gs[1, 0:3])


fig = gcf()
gs.tight_layout(fig,h_pad=12,w_pad=1)
ax = [ax1, ax2]



ax[0].bar(x - width*2, enso_control, width,label='PiControl', color='slategray')
ax[0].bar(x - width, enso_ssp585_0, width,label='200yr',color='royalblue')
ax[0].bar(x, enso_ssp585_1, width,label='400yr',color='mediumpurple')
ax[0].bar(x + width, enso_ssp585_2, width,label='600yr',color='darkorchid')
ax[0].bar(x + width*2, enso_ssp585_3, width,label='800yr',color='blue')
ax[0].bar(x + width*3, enso_ssp585_4, width,label='1000yr',color='navy')

ax[0].set_xticks(range(len(model_names)))

ax[0].set_xticklabels(model_names,rotation='horizontal')

ax[0].set_title('ENSO in CMIP6 long simulations',fontsize=40)

ax[0].set_ylabel('ENSO SST amp. ($^o$C)')
ax[0].legend(ncol=2,fontsize=38)

ax[1].bar(x - width*2, noise_control, width,label='PiControl', color='brown')
ax[1].bar(x - width, noise_ssp585_0, width,label='200yr', color='pink')
ax[1].bar(x, noise_ssp585_1, width,label='400yr', color='grey')
ax[1].bar(x + width, noise_ssp585_2, width,label='600yr', color='olive')
ax[1].bar(x + width*2, noise_ssp585_3, width,label='800yr', color='cyan')
ax[1].bar(x + width*3, noise_ssp585_4, width,label='1000yr',color='orange')

ax[1].set_xticks(range(len(model_names)))
ax[1].set_xticklabels(model_names,rotation='vertical')

ax[1].set_title('noise in CMIP6 long simulations',fontsize=40)

ax[1].set_ylabel('eq. std')
ax[1].legend(ncol=2,fontsize=38)

#plt.figtext(0.01, 0.97, 'a)',fontweight='normal')
#plt.figtext(0.4, 0.97, 'b)',fontweight='normal')
#plt.figtext(0.6, 0.97, 'c)',fontweight='normal')
#%%
