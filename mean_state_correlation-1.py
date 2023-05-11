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

model_names=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5-CanOE','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']

model_names_1=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']
    
# model_names_1=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2-WACCM','CESM2',\
#              'EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
#                  'HadGEM3-GC31-LL','HadGEM3-GC31-MM','KACE-1-0-G',\
#                     'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']
model_names_2=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'IPSL-CM6A-LR','KACE-1-0-G',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']
    
#model_names=['CanESM5','CESM2-WACCM','CESM2',\
#             'CNRM-CM6-1','CNRM-ESM2-1','FGOALS-g3',\
#                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
#                     'MPI-ESM1-2-LR','UKESM1-0-LL']



nino3_start=180
nino3_stop=280

nino4_start=160
nino4_stop=210
    
e1=180
e2=280

w1=150
w2=180
x=20

grad_dif_4x=[0]*len(model_names_1)
grad_dif_1pct=[0]*len(model_names)
grad_dif_ssp5=[0]*len(model_names)
grad_dif_ssp1=[0]*len(model_names_2)

enso_4x=[0]*len(model_names_1)
enso_control=[0]*len(model_names)
enso_ssp585_1=[0]*len(model_names)
enso_ssp585_2=[0]*len(model_names)
enso_ssp585_3=[0]*len(model_names)
ssp_total=[0]*len(model_names)
ssp_max=[0]*len(model_names)
ssp_min=[0]*len(model_names)
ssp_total_126=[0]*len(model_names_2)
ssp_max_126=[0]*len(model_names_2)
ssp_min_126=[0]*len(model_names_2)
enso_ssp126=[0]*len(model_names_2)
enso_1pct=[0]*len(model_names)

#%%
for x in range(0,len(model_names_1)):
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq4_'+model_names_1[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(50*12,150*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')

    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_1[x]+'.nc')
    control=control['ts']
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    
    grad_dif_4x[x]=(grad1-gradC).values

    enso_ssp=(enso1)
    enso_4x[x]=enso1
    
for x in range(0,len(model_names)):
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names[x]+'.nc')
    
    control=control['ts']
    
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    control=control.groupby('time.month')-control.groupby('time.month').mean('time')
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=control.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    ensoC=low_nino3.std('time')
    
    
    enso_control[x]=ensoC
    
  
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1pct_'+model_names[x]+'.nc')
    
    control=control['ts'].isel(time=slice(50*12,150*12))
    
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    control=control.groupby('time.month')-control.groupby('time.month').mean('time')
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=control.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')
    
    
    enso_1pct[x]=enso1
    
    grad_dif_1pct[x]=(grad1-gradC).values
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1_'+model_names[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_1[x]=enso1
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq2_'+model_names[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad2=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso2=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_2[x]=enso2
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq3_'+model_names[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad3=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso3=low_nino3.std('time')


    grad_ssp=(grad1+grad2+grad3)/3
    grad_dif_ssp5[x]=(grad_ssp-gradC).values
    
    enso_ssp=(enso1)
    enso_ssp585_3[x]=enso3
    test=xr.concat([enso1,enso2,enso3], dim='ens')
    ssp_total[x]=test.mean('ens')
    ssp_max[x]=test.max('ens')-ssp_total[x]
    ssp_min[x]=ssp_total[x]-test.min('ens')


    plt.plot(ensoC)
    plt.plot(enso_ssp)
#plt.plot(control.mean('lat').mean('time'))
#plt.plot(ens1.mean('lat').mean('time'))

#%%
for x in range(0,len(model_names_2)):
    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_2[x]+'.nc')
    
    control=control['ts']
    
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    control=control.groupby('time.month')-control.groupby('time.month').mean('time')
    control=control.assign_coords(time=list(range(len(control.time))))
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1ssp126_'+model_names_2[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_1[x]=enso1
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq2ssp126_'+model_names_2[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad3=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso2=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_2[x]=enso2
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq3ssp126_'+model_names_2[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,85*12))

    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad3=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso3=low_nino3.std('time')


    grad_ssp=(grad1+grad2+grad3)/3
    grad_dif_ssp1[x]=(grad_ssp-gradC).values
    enso_ssp=(enso1)
    enso_ssp585_3[x]=enso3
    test=xr.concat([enso1,enso2,enso3], dim='ens')
    ssp_total_126[x]=test.mean('ens')
    ssp_max_126[x]=test.max('ens')-ssp_total_126[x]
    ssp_min_126[x]=ssp_total_126[x]-test.min('ens')
#%%
enso_4x_pr=[0]*len(model_names_1)
enso_control_pr=[0]*len(model_names)
enso_ssp585_1=[0]*len(model_names)
enso_ssp585_2=[0]*len(model_names)
enso_ssp585_3=[0]*len(model_names)
ssp_total_pr=[0]*len(model_names)
ssp_max_pr=[0]*len(model_names)
ssp_min_pr=[0]*len(model_names)
ssp_total_126_pr=[0]*len(model_names_2)
ssp_max_126_pr=[0]*len(model_names_2)
ssp_min_126_pr=[0]*len(model_names_2)
enso_ssp126=[0]*len(model_names_2)
enso_1pct_pr=[0]*len(model_names)

#%%
for x in range(0,len(model_names_1)):
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/pr_eq4_'+model_names_1[x]+'.nc')
    
    ens1=ens1['pr'].isel(time=slice(50*12,150*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_4x_pr[x]=enso1*86400
    
for x in range(0,len(model_names)):
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/pr_eqC_'+model_names[x]+'.nc')
    control=control['pr']
    
    
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    control=control.groupby('time.month')-control.groupby('time.month').mean('time')
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=control.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    ensoC=low_nino3.std('time')
    
    
    enso_control_pr[x]=ensoC*86400
    
  
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/pr_eqpct_'+model_names[x]+'.nc')
    control=control['pr'].isel(time=slice(50*12,150*12))
 
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    control=control.groupby('time.month')-control.groupby('time.month').mean('time')
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=control.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')
    
    
    enso_1pct_pr[x]=enso1*86400
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/pr_eq1_'+model_names[x]+'.nc')
    ens1=ens1['pr'].isel(time=slice(0,85*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))   
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_1[x]=enso1
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/pr_eq2_'+model_names[x]+'.nc')
    ens1=ens1['pr'].isel(time=slice(0,85*12))

    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso2=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_2[x]=enso2
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/pr_eq3_'+model_names[x]+'.nc')
    ens1=ens1['pr'].isel(time=slice(0,85*12))

    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso3=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_3[x]=enso3
    test=xr.concat([enso1,enso2,enso3], dim='ens')*86400
    ssp_total_pr[x]=test.mean('ens')
    ssp_max_pr[x]=(test.max('ens')-ssp_total_pr[x])
    ssp_min_pr[x]=(ssp_total_pr[x]-test.min('ens'))


#plt.plot(control.mean('lat').mean('time'))
#plt.plot(ens1.mean('lat').mean('time'))

#%%
for x in range(0,len(model_names_2)):
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/pr_eq1ssp126_'+model_names_2[x]+'.nc')
    ens1=ens1['pr'].isel(time=slice(0,85*12))

    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_1[x]=enso1
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/pr_eq2ssp126_'+model_names_2[x]+'.nc')
    ens1=ens1['pr'].isel(time=slice(0,85*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso2=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_2[x]=enso2
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/pr_eq3ssp126_'+model_names_2[x]+'.nc')
    ens1=ens1['pr'].isel(time=slice(0,85*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso3=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_3[x]=enso3
    test=xr.concat([enso1,enso2,enso3], dim='ens')*86400
    ssp_total_126_pr[x]=test.mean('ens')
    ssp_max_126_pr[x]=(test.max('ens')-ssp_total_126_pr[x])
    ssp_min_126_pr[x]=(ssp_total_126_pr[x]-test.min('ens'))

#%%
plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'hatch.color': '0.1'})  

width = 0.15  # the width of the bars
fig = figure(figsize=(47,12.5))
gs = gridspec.GridSpec(1, 3)
ax1 = plt.subplot(gs[0, 0:3])
#ax2 = plt.subplot(gs[1, 0:1])
#ax3 = plt.subplot(gs[1, 1:2])
#ax4 = plt.subplot(gs[1, 2:3])

fig1 = gcf()
gs.tight_layout(fig1,h_pad=12,w_pad=1)
ax = [ax1]
x = np.arange(len(model_names))  # the label locations
x1=array([ 0,  1, 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,16,17,18])
#x1=array([ 0,  1, 3,  4,    7,  8,  9,  11, 12, 13, 14, 15,16,17,18])
x2=array([ 0,  1,  3,  5,  6,  7,  8,  9, 10,11, 14, 15,16,17,18])


ax[0].bar(x - width*2, enso_control_pr, width,label='PiControl',color='k')
ax[0].bar(x1 - width, enso_4x_pr, width,label='abrupt-4xCO2',color='orange')
ax[0].bar(x, enso_1pct_pr, width,label='1pctCO2',color='green')
ax[0].bar(x + width, ssp_total_pr, width,yerr=[ssp_min_pr,ssp_max_pr], error_kw=dict(lw=5),label='ssp585',color='red')
ax[0].bar(x2 + width*2, ssp_total_126_pr, width,yerr=[ssp_min_126_pr,ssp_max_126_pr], error_kw=dict(lw=5),label='ssp126',color='purple')

ax[0].bar(20 - width*2,mean(enso_control_pr),width, color='k')
ax[0].bar(20 - width,mean(enso_4x_pr),width, color='orange')
ax[0].bar(20,mean(enso_1pct_pr),width, color='green')
ax[0].bar(20 + width,mean(ssp_total_pr),width, color='red')
ax[0].bar(20 + width*2,mean(ssp_total_126_pr),width, color='purple')

model_names3=model_names+['']+['multi-model mean']

ax[0].set_xticks(range(len(model_names3)))

ax[0].set_xticklabels(model_names3,rotation='vertical')


ax[0].set_xticks(range(len(model_names3)))


ax[0].set_title('ENSO rainfall in CMIP6 models',fontsize=40)

ax[0].set_ylabel('ENSO rainfall amp. (mm/day)')
ax[0].legend(ncol=2,fontsize=38)
#%%

enso_control_pr=np.array(enso_control_pr)

enso_control_1_pr=enso_control_pr[x1]
enso_control_2_pr=enso_control_pr[x2]

enso_4x_pr=np.array(enso_4x_pr)
enso_1pct_pr=np.array(enso_1pct_pr)


ssp_total_pr=np.array(ssp_total_pr)
ssp_total_126_pr=np.array(ssp_total_126_pr)

diff_4x_pr = enso_4x_pr-enso_control_1_pr
 
diff_1pct_pr=enso_1pct_pr-enso_control_pr
diff_ssp5_pr=ssp_total_pr-enso_control_pr
diff_ssp1_pr=ssp_total_126_pr-enso_control_2_pr

markerlist=np.array(['o','v','^','<','>','1','2','3','4','s','p','*','h','+','x','X','D','|','_','o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_'])
colorlist=np.array(['grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k','brown','orange','olive','green','cyan','blue','purple','pink','red','k',\
    'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k', 'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k'])


#%%
plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'hatch.color': '0.1'})  
x = np.arange(len(model_names))  # the label locations
x1=array([ 0,  1, 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,16,17,18])
#x1=array([ 0,  1, 3,  4,    7,  8,  9,  11, 12, 13, 14, 15,16,17,18])
x2=array([ 0,  1,  3,  5,  6,  7,  8,  9, 10,11, 14, 15,16,17,18])


width = 0.15  # the width of the bars
fig2 = figure(figsize=(47,28))
gs = gridspec.GridSpec(2, 3)

ax2 = plt.subplot(gs[0, 0:1])
ax3 = plt.subplot(gs[0, 1:2])
ax4 = plt.subplot(gs[0, 2:3])
ax5 = plt.subplot(gs[1, 0:1])
ax6 = plt.subplot(gs[1, 1:2])
ax7 = plt.subplot(gs[1, 2:3])
#ax8 = plt.subplot(gs[2, 0:1])
#ax9 = plt.subplot(gs[2, 1:2])
#ax10 = plt.subplot(gs[2, 2:3])
fig2 = gcf()
gs.tight_layout(fig2,h_pad=4,w_pad=4)
ax = [ax2, ax3, ax4, ax5, ax6, ax7]

plt.figtext(0.04, 0.975, 'a)')
plt.figtext(0.372, 0.975, 'b)')
plt.figtext(0.72, 0.975, 'c)')

plt.figtext(0.04, 0.46, 'd)')
plt.figtext(0.38, 0.46, 'e)')

plt.figtext(0.72, 0.46, 'f)')

enso_control=np.array(enso_control)

enso_control_1=enso_control[x1]
enso_control_2=enso_control[x2]

enso_4x=np.array(enso_4x)
enso_1pct=np.array(enso_1pct)


ssp_total=np.array(ssp_total)
ssp_total_126=np.array(ssp_total_126)

diff_4x = enso_4x-enso_control_1
 
diff_1pct=enso_1pct-enso_control
diff_ssp5=ssp_total-enso_control
diff_ssp1=ssp_total_126-enso_control_2

diff_4x_p =( enso_4x-enso_control_1)/enso_control_1*100
 
diff_1pct_p=(enso_1pct-enso_control)/enso_control*100
diff_ssp5_p=(ssp_total-enso_control)/enso_control*100
diff_ssp1_p=(ssp_total_126-enso_control_2)/enso_control_2*100
#fuk this fucking shitturd

markerlist_subset1=markerlist[x1]

colorlist_subset1=colorlist[x1]

for i in x:
    ax[0].scatter(grad_dif_1pct[i],diff_1pct_pr[i],s=1000,marker=markerlist[i],c=colorlist[i])

for i in range(0,len(x1)):
    ax[1].scatter(grad_dif_4x[i],diff_4x_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])
    ax[2].scatter(grad_dif_ssp5[i],diff_ssp5_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])


ax[0].set_title('1pctCO$_2$ mean state vs ENSO rainfall',fontsize=40)
ax[0].set_xlabel('$\Delta$gradient ($^o$C)')
ax[0].set_ylabel('$\Delta$std rainfall (mm/day)')
ax[0].set_ylim(0.2*10**(-0),1*10**(-0))

ax[1].set_title('abrupt4xCO$_2$ mean state vs ENSO rainf.',fontsize=40)
ax[1].set_xlabel('$\Delta$gradient ($^o$C)')
ax[1].set_ylabel('$\Delta$std rainfall (mm/day)')
ax[1].set_ylim(0.2*10**(-0),2*10**(-0))

ax[2].set_title('SSP585 mean state vs ENSO rainfall',fontsize=40)
ax[2].set_xlabel('$\Delta$gradient ($^o$C)')
ax[2].set_ylabel('$\Delta$std rainfall (mm/day)')
ax[2].set_ylim(0.2*10**(-0),2*10**(-0))

# for i in x:
#     ax[3].scatter(grad_dif_1pct[i],diff_1pct[i],s=800,marker=markerlist[i],c=colorlist[i])
# for i in range(0,len(x1)):
#     ax[4].scatter(grad_dif_4x[i],diff_4x[i],s=800,marker=markerlist_subset1[i],c=colorlist_subset1[i])
#     ax[5].scatter(grad_dif_ssp5[i],diff_ssp5[i],s=800,marker=markerlist_subset1[i],c=colorlist_subset1[i])


# ax[3].set_title('1pctCO2 mean state vs ENSO SST',fontsize=40)
# ax[3].set_xlabel('$\Delta$gradient')
# ax[3].set_ylabel('$\Delta$std SST')
# ax[3].set_ylim(-0.2,0.55)

# ax[4].set_title('4xCO2 mean state vs ENSO SST',fontsize=40)
# ax[4].set_xlabel('$\Delta$gradient')
# ax[4].set_ylabel('$\Delta$std SST')
# ax[4].set_ylim(-0.2,0.55)

# ax[5].set_title('ssp585 mean state vs ENSO SST',fontsize=40)
# ax[5].set_xlabel('$\Delta$gradient')
# ax[5].set_ylabel('$\Delta$std SST')
# ax[5].set_ylim(-0.2,0.55)

for i in x:
    ax[3].scatter(diff_1pct[i],diff_1pct_pr[i],s=1000,marker=markerlist[i],c=colorlist[i])
for i in range(0,len(x1)):
    ax[4].scatter(diff_4x[i],diff_4x_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])
    ax[5].scatter(diff_ssp5[i],diff_ssp5_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])

ax[3].set_title('1pctCO$_2$ ENSO SST vs ENSO rainfall',fontsize=40)
ax[3].set_xlabel('$\Delta$std SST ($^o$C)')
ax[3].set_ylabel('$\Delta$std rainfall (mm/day)')
ax[3].set_ylim(0.2*10**(-0),1*10**(-0))

ax[4].set_title('4xCO$_2$ ENSO SST vs ENSO rainfall',fontsize=40)
ax[4].set_xlabel('$\Delta$std SST ($^o$C)')
ax[4].set_ylabel('$\Delta$std rainfall (mm/day)')
ax[4].set_ylim(0.2*10**(-0),2*10**(-0))

ax[5].set_title('SSP585 ENSO SST vs ENSO rainfall',fontsize=40)
ax[5].set_xlabel('$\Delta$std SST ($^o$C)')
ax[5].set_ylabel('$\Delta$std rainfall (mm/day)')
ax[5].set_ylim(0.2*10**(-0),2*10**(-0))

from scipy.stats import kendalltau, pearsonr, spearmanr

def pearsonr_pval(x,y):
        return pearsonr(x,y)[1]
def pearsonr_r(x,y):
        return pearsonr(x,y)[0]

diff_1pctS=pd.Series(diff_1pct)
diff_ssp5S=pd.Series(diff_ssp5)
diff_4xS=pd.Series(diff_4x)

diff_1pct_prS=pd.Series(diff_1pct_pr)
diff_ssp5_prS=pd.Series(diff_ssp5_pr)
diff_4x_prS=pd.Series(diff_4x_pr)

grad_dif_1pctS=pd.Series(grad_dif_1pct)
grad_dif_ssp5S=pd.Series(grad_dif_ssp5)
grad_dif_4xS=pd.Series(grad_dif_4x)

corr1=grad_dif_1pctS.corr(diff_1pct_prS,method=pearsonr_r)
pval1 = grad_dif_1pctS.corr(diff_1pct_prS,method=pearsonr_pval)

corr2=grad_dif_4xS.corr(diff_4x_prS,method=pearsonr_r)
pval2 = grad_dif_4xS.corr(diff_4x_prS,method=pearsonr_pval)

corr3=grad_dif_ssp5S.corr(diff_ssp5_prS,method=pearsonr_r)
pval3 = grad_dif_ssp5S.corr(diff_ssp5_prS,method=pearsonr_pval)

corr4=grad_dif_1pctS.corr(diff_1pctS,method=pearsonr_r)
pval4 = grad_dif_1pctS.corr(diff_1pctS,method=pearsonr_pval)

corr5=grad_dif_4xS.corr(diff_4xS,method=pearsonr_r)
pval5 = grad_dif_4xS.corr(diff_4xS,method=pearsonr_pval)

corr6=grad_dif_ssp5S.corr(diff_ssp5S,method=pearsonr_r)
pval6 = grad_dif_ssp5S.corr(diff_ssp5S,method=pearsonr_pval)

corr7=diff_1pctS.corr(diff_1pct_prS,method=pearsonr_r)
pval7 = diff_1pctS.corr(diff_1pctS,method=pearsonr_pval)

corr8=diff_4xS.corr(diff_4x_prS,method=pearsonr_r)
pval8 = diff_4xS.corr(diff_4xS,method=pearsonr_pval)

corr9=diff_ssp5S.corr(diff_ssp5_prS,method=pearsonr_r)
pval9 = diff_ssp5S.corr(diff_ssp5_prS,method=pearsonr_pval)

plt.figtext(0.05,0.92,'R='+str("%.2f" % corr1)+',  p='+str("%.2f" % pval1))
plt.figtext(0.4,0.92,'R='+str("%.2f" % corr2)+',  p='+str("%.2f" % pval2))
plt.figtext(0.85,0.92,'R='+str("%.2f" % corr3)+',  p='+str("%.2f" % pval3))

#plt.figtext(0.05,0.37,'R='+str("%.2f" % corr4)+',  p='+str("%.2f" % pval4))
#plt.figtext(0.4,0.37,'R='+str("%.2f" % corr5)+',  p='+str("%.2f" % pval5))
#plt.figtext(0.85,0.37,'R='+str("%.2f" % corr6)+',  p='+str("%.2f" % pval6))

plt.figtext(0.05,0.40,'R='+str("%.2f" % corr7)+',  p='+str("%.2f" % pval7))
plt.figtext(0.4,0.40,'R='+str("%.2f" % corr8)+',  p='+str("%.2f" % pval8))
plt.figtext(0.75,0.40,'R='+str("%.2f" % corr9)+',  p='+str("%.2f" % pval9))
#%%



#%%
plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'hatch.color': '0.1'})  
x = np.arange(len(model_names))  # the label locations
x1=array([ 0,  1, 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,16,17,18])
x2=array([ 0,  1,  3,  5,  6,  7,  8,  9, 10,11, 14, 15,16,17,18])


width = 0.15  # the width of the bars
fig3 = figure(figsize=(47,26))
gs = gridspec.GridSpec(2, 3)

ax2 = plt.subplot(gs[0, 0:1])
ax3 = plt.subplot(gs[0, 1:2])
ax4 = plt.subplot(gs[0, 2:3])

ax5 = plt.subplot(gs[1, 0:1])
ax6 = plt.subplot(gs[1, 1:2])
ax7 = plt.subplot(gs[1, 2:3])

#ax8 = plt.subplot(gs[2, 0:1])
#ax9 = plt.subplot(gs[2, 1:2])
#ax10 = plt.subplot(gs[2, 2:3])
fig23= gcf()
gs.tight_layout(fig3,h_pad=4,w_pad=4)
ax = [ax2, ax3, ax4, ax5, ax6, ax7]



# for i in x:
#     ax[0].scatter(grad_dif_1pct[i],diff_1pct_pr[i],s=1000,marker=markerlist[i],c=colorlist[i])

# for i in range(0,len(x1)):
#     ax[1].scatter(grad_dif_4x[i],diff_4x_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])
#     ax[2].scatter(grad_dif_ssp5[i],diff_ssp5_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])


# ax[0].set_title('1pctCO2 mean state vs ENSO rainfall',fontsize=40)
# ax[0].set_xlabel('$\Delta$gradient')
# ax[0].set_ylabel('$\Delta$std rainfall')
# ax[0].set_ylim(0.2*10**(-5),2*10**(-5))

# ax[1].set_title('abrupt4xCO2 mean state vs ENSO rainfall',fontsize=40)
# ax[1].set_xlabel('$\Delta$gradient')
# ax[1].set_ylabel('$\Delta$std rainfall')
# ax[1].set_ylim(0.2*10**(-5),2*10**(-5))

# ax[2].set_title('SSP585 mean state vs ENSO rainfall',fontsize=40)
# ax[2].set_xlabel('$\Delta$gradient')
# ax[2].set_ylabel('$\Delta$std rainfall')
# ax[2].set_ylim(0.2*10**(-5),2*10**(-5))

for i in x:
    ax[0].scatter(grad_dif_1pct[i],diff_1pct[i],s=1000,marker=markerlist[i],c=colorlist[i])
for i in range(0,len(x1)):
    ax[1].scatter(grad_dif_4x[i],diff_4x[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])
    ax[2].scatter(grad_dif_ssp5[i],diff_ssp5[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])


ax[0].set_title('1pctCO$_2$: ENSO amp. vs mean grad.',fontsize=40)
ax[0].set_xlabel('$\Delta$gradient ($^o$C)')
ax[0].set_ylabel('$\Delta$std SST ($^o$C)')
ax[0].set_ylim(-0.2,0.6)

ax[1].set_title('abrupt4xCO$_2$: ENSO amp. vs mean grad.',fontsize=40)
ax[1].set_xlabel('$\Delta$gradient ($^o$C)')
ax[1].set_ylabel('$\Delta$std SST ($^o$C)')
ax[1].set_ylim(-0.2,0.6)

ax[2].set_title('ssp585: ENSO amp. vs mean grad.',fontsize=40)
ax[2].set_xlabel('$\Delta$gradient ($^o$C)')
ax[2].set_ylabel('$\Delta$std SST ($^o$C)')
ax[2].set_ylim(-0.2,0.6)

# for i in x:
#     ax[3].scatter(diff_1pct[i],diff_1pct_pr[i],s=1000,marker=markerlist[i],c=colorlist[i])
# for i in range(0,len(x1)):
#     ax[4].scatter(diff_4x[i],diff_4x_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])
#     ax[5].scatter(diff_ssp5[i],diff_ssp5_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])

# ax[3].set_title('1pctCO2 ENSO SST vs ENSO rainfall',fontsize=40)
# ax[3].set_xlabel('$\Delta$std SST')
# ax[3].set_ylabel('$\Delta$std rainfall')
# ax[3].set_ylim(0.2*10**(-5),2*10**(-5))

# ax[4].set_title('4xCO2 ENSO SST vs ENSO rainfall',fontsize=40)
# ax[4].set_xlabel('$\Delta$std SST')
# ax[4].set_ylabel('$\Delta$std rainfall')
# ax[4].set_ylim(0.2*10**(-5),2*10**(-5))

# ax[5].set_title('SSP585 ENSO SST vs ENSO rainfall',fontsize=40)
# ax[5].set_xlabel('$\Delta$std SST')
# ax[5].set_ylabel('$\Delta$std rainfall')
# ax[5].set_ylim(0.2*10**(-5),2*10**(-5))



# plt.figtext(0.05,0.92,'R='+str("%.2f" % corr1)+',  p='+str("%.2f" % pval1))
# plt.figtext(0.4,0.92,'R='+str("%.2f" % corr2)+',  p='+str("%.2f" % pval2))
# plt.figtext(0.85,0.92,'R='+str("%.2f" % corr3)+',  p='+str("%.2f" % pval3))


# plt.figtext(0.05,0.40,'R='+str("%.2f" % corr7)+',  p='+str("%.2f" % pval7))
# plt.figtext(0.4,0.40,'R='+str("%.2f" % corr8)+',  p='+str("%.2f" % pval8))
# plt.figtext(0.75,0.40,'R='+str("%.2f" % corr9)+',  p='+str("%.2f" % pval9))



#%%
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

model_names=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5-CanOE','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']

x=([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])


x1=array([0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,18,19])

x2=array([0,1,3,4,5,6,7,12,14,16,17,18,19])

x3=array([0,1,3,4,5,6,7,12,16,17,18,19])

model_names_4x=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1',\
                 'HadGEM3-GC31-LL','IPSL-CM6A-LR',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']
    
model_names_ssp=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5-CanOE','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','HadGEM3-GC31-LL','IPSL-CM6A-LR',\
                    'MIROC6', 'MPI-ESM1-2-LR','UKESM1-0-LL']


   
model_names_1pct=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1',\
                 'HadGEM3-GC31-LL',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']




nino3_start=180
nino3_stop=280

nino4_start=160
nino4_stop=210
    
e1=180
e2=280

w1=150
w2=180
x=20

grad_dif_4x=[0]*len(model_names_4x)
grad_dif_1pct=[0]*len(model_names_1pct)
grad_dif_ssp5=[0]*len(model_names_ssp)
#grad_dif_ssp1=[0]*len(model_names_2)

strat_4x=[0]*len(model_names_4x)
strat_1pct=[0]*len(model_names_1pct)
strat_ssp5_1=[0]*len(model_names_ssp)
strat_ssp5_2=[0]*len(model_names_ssp)
strat_ssp5_3=[0]*len(model_names_ssp)

#strat_ssp1=[0]*len(model_names_2)

enso_4x=[0]*len(model_names_4x)
enso_1pct=[0]*len(model_names_1pct)
enso_control=[0]*len(model_names_ssp)
enso_ssp585_1=[0]*len(model_names_ssp)
enso_ssp585_2=[0]*len(model_names_ssp)
enso_ssp585_3=[0]*len(model_names_ssp)
ssp_total=[0]*len(model_names_ssp)
ssp_max=[0]*len(model_names_ssp)
ssp_min=[0]*len(model_names_ssp)

ssp_total_strat=[0]*len(model_names_ssp)
ssp_max_strat=[0]*len(model_names_ssp)
ssp_min_strat=[0]*len(model_names_ssp)

strat_4x_lon = np.zeros((71,len(model_names_4x)))
strat_1pct_lon = np.zeros((71,len(model_names_1pct)))
strat_ssp_lon_1 = np.zeros((71,len(model_names_ssp)))
strat_ssp_lon_2 = np.zeros((71,len(model_names_ssp)))
strat_ssp_lon_3 = np.zeros((71,len(model_names_ssp)))
strat_ssp_lon_mean = np.zeros((71,len(model_names_ssp)))

#%%
for x in range(0,len(model_names_4x)):
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq4_'+model_names_4x[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,150*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')

    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_4x[x]+'.nc')
    control=control['ts']
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    ens1=control.groupby('time.month')-control.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    ensoC=low_nino3.std('time')
    
    grad_dif_4x[x]=(grad1-gradC).values

    enso_ssp=(enso1)
    enso_4x[x]=(enso1-ensoC).values
    
    ens1_strat=xr.open_dataset('/Users/ullaheede_1/Downloads/strat_'+model_names_4x[x]+'_4x.nc')
    pi_strat=xr.open_dataset('/Users/ullaheede_1/Downloads/strat_'+model_names_4x[x]+'_piControl.nc')
    strat_4x[x]=(ens1_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)).values
    strat_4x_lon[:,x]=(ens1_strat['thetao'].mean('lat',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True)).squeeze().values
    
for x in range(0,len(model_names_1pct)):
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1pct_'+model_names_1pct[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,150*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')

    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_1pct[x]+'.nc')
    control=control['ts']
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    ens1=control.groupby('time.month')-control.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    ensoC=low_nino3.std('time')
    
    grad_dif_1pct[x]=(grad1-gradC).values

    enso_ssp=(enso1)
    enso_1pct[x]=(enso1-ensoC).values
    
    ens1_strat=xr.open_dataset('/Users/ullaheede_1/Downloads/strat_'+model_names_1pct[x]+'_1pct.nc')
    pi_strat=xr.open_dataset('/Users/ullaheede_1/Downloads/strat_'+model_names_1pct[x]+'_piControl.nc')
    strat_1pct[x]=(ens1_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)).values
    strat_1pct_lon[:,x]=(ens1_strat['thetao'].mean('lat',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True)).squeeze().values
    
for x in range(0,len(model_names_ssp)):
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_ssp[x]+'.nc')
    
    control=control['ts']
    
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    control=control.groupby('time.month')-control.groupby('time.month').mean('time')
    control=control.assign_coords(time=list(range(len(control.time))))
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=control.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    ensoC=low_nino3.std('time')
    pi_strat=xr.open_dataset('/Users/ullaheede_1/Downloads/strat_'+model_names_ssp[x]+'_piControl.nc')
     
    
    enso_control[x]=ensoC
    
  
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1_'+model_names_ssp[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')
    ens1_strat=xr.open_dataset('/Users/ullaheede_1/Downloads/strat_'+model_names_ssp[x]+'_ssp585_1.nc')
    strat_ssp5_1[x]=ens1_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)
    strat_ssp_lon_1[:,x]=(ens1_strat['thetao'].mean('lat',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True)).squeeze().values
 
    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_1[x]=enso1
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq2_'+model_names_ssp[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad2=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso2=low_nino3.std('time')
    ens1_strat=xr.open_dataset('/Users/ullaheede_1/Downloads/strat_'+model_names_ssp[x]+'_ssp585_2.nc')
    strat_ssp5_2[x]=ens1_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)
    strat_ssp_lon_2[:,x]=(ens1_strat['thetao'].mean('lat',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True)).squeeze().values


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_2[x]=enso2
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq3_'+model_names_ssp[x]+'.nc')
    
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    
    
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad3=(west-east).mean('time')
    ens1=ens1#.groupby('time.month')-ens1.groupby('time.month').mean('time')
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso3=low_nino3.std('time')
    ens1_strat=xr.open_dataset('/Users/ullaheede_1/Downloads/strat_'+model_names_ssp[x]+'_ssp585_3.nc')
    strat_ssp5_3[x]=ens1_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True).mean('lon',skipna=True)
    strat_ssp_lon_3[:,x]=(ens1_strat['thetao'].mean('lat',skipna=True)-pi_strat['thetao'].mean('lat',skipna=True)).squeeze().values


    grad_ssp=(grad1+grad2+grad3)/3
    grad_dif_ssp5[x]=(grad_ssp-gradC).values
    
    enso_ssp=(enso1)
    enso_ssp585_3[x]=enso3
    test=xr.concat([enso1,enso2,enso3], dim='ens')
    ssp_total[x]=test.mean('ens')
    ssp_max[x]=test.max('ens')-ssp_total[x]
    ssp_min[x]=ssp_total[x]-test.min('ens')
    
    test=xr.concat([strat_ssp5_1[x],strat_ssp5_2[x],strat_ssp5_3[x]], dim='ens')
    ssp_total_strat[x]=(test.mean('ens')).values
    ssp_max_strat[x]=test.max('ens')-ssp_total_strat[x]
    ssp_min_strat[x]=ssp_total_strat[x]-test.min('ens')
    strat_ssp_lon_mean[:,x]=(strat_ssp_lon_1[:,x]+strat_ssp_lon_2[:,x]+strat_ssp_lon_3[:,x])/3

    #test=xr.concat([enso1,enso2,enso3], dim='ens')
    #ssp_total[x]=test.mean('ens')
    #ssp_max[x]=test.max('ens')-ssp_total[x]
    #ssp_min[x]=ssp_total[x]-test.min('ens')


    plt.plot(ensoC)
    plt.plot(enso_ssp)
#plt.plot(control.mean('lat').mean('time'))
#plt.plot(ens1.mean('lat').mean('time'))
ssp_total=np.array(ssp_total)
enso_control=np.array(enso_control)
diff_ssp5=ssp_total-enso_control

plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'hatch.color': '0.1'})  
x = np.arange(len(model_names))  # the label locations

markerlist=np.array(['o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_','o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_'])
colorlist=np.array(['grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k','grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k',\
    'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k', 'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k'])

x=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

x1=array([0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,18,19])

x2=array([0,1,3,4,5,6,7,12,14,16,17,18,19])

x3=array([0,1,3,4,5,6,7,12,16,17,18,19])

markerlist_subset1=markerlist[x1]
markerlist_subset2=markerlist[x2]
markerlist_subset3=markerlist[x3]
colorlist_subset1=colorlist[x1]
colorlist_subset2=colorlist[x2]
colorlist_subset3=colorlist[x3]

width = 0.15  # the width of the bars


plt.figtext(0.04, 0.975, 'a)')
plt.figtext(0.372, 0.975, 'b)')
plt.figtext(0.72, 0.975, 'c)')
plt.figtext(0.04, 0.46, 'e)')
plt.figtext(0.372, 0.46, 'f)')
plt.figtext(0.72, 0.46, 'g)')
# for i in x:
#     ax[0].scatter(grad_dif_1pct[i],diff_1pct_pr[i],s=1000,marker=markerlist[i],c=colorlist[i])

# for i in range(0,len(x1)):
#     ax[1].scatter(grad_dif_4x[i],diff_4x_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])
#     ax[2].scatter(grad_dif_ssp5[i],diff_ssp5_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])


# ax[0].set_title('1pctCO2 mean state vs ENSO rainfall',fontsize=40)
# ax[0].set_xlabel('$\Delta$gradient')
# ax[0].set_ylabel('$\Delta$std rainfall')
# ax[0].set_ylim(0.2*10**(-5),2*10**(-5))

# ax[1].set_title('abrupt4xCO2 mean state vs ENSO rainfall',fontsize=40)
# ax[1].set_xlabel('$\Delta$gradient')
# ax[1].set_ylabel('$\Delta$std rainfall')
# ax[1].set_ylim(0.2*10**(-5),2*10**(-5))

# ax[2].set_title('SSP585 mean state vs ENSO rainfall',fontsize=40)
# ax[2].set_xlabel('$\Delta$gradient')
# ax[2].set_ylabel('$\Delta$std rainfall')
# ax[2].set_ylim(0.2*10**(-5),2*10**(-5))

#for i in x:
#    ax[0].scatter(grad_dif_1pct[i],diff_1pct[i],s=1000,marker=markerlist[i],c=colorlist[i])
for i in range(0,len(x2)):
    ax[4].scatter(strat_4x[i],enso_4x[i],s=1000,marker=markerlist_subset2[i],c=colorlist_subset2[i])

for i in range(0,len(x1)):
    ax[5].scatter(ssp_total_strat[i],diff_ssp5[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])
for i in range(0,len(x3)):
    ax[3].scatter(strat_1pct[i],enso_1pct[i],s=1000,marker=markerlist_subset3[i],c=colorlist_subset3[i])
    #  ax[3].scatter(grad_dif_ssp5[i],ssp_total_strat[i],s=1000,marker=markerlist[i],c=colorlist[i])
  #  ax[3].scatter(grad_dif_4x[i],strat_4x[i],s=1000,marker=markerlist[i],c=colorlist[i])
  #  ax[3].scatter(grad_dif_1pct[i],strat_1pct[i],s=1000,marker=markerlist[i],c=colorlist[i])

strat_4xS=pd.Series(strat_4x)
ssp_total_stratS=pd.Series(ssp_total_strat)
strat_1pctS=pd.Series(strat_1pct)
enso_4xS=pd.Series(enso_4x)
diff_ssp5S=pd.Series(diff_ssp5)
enso_1pctS=pd.Series(enso_1pct)

from scipy.stats import kendalltau, pearsonr, spearmanr

def pearsonr_pval(x,y):
        return pearsonr(x,y)[1]
def pearsonr_r(x,y):
        return pearsonr(x,y)[0]

corr1=strat_1pctS.corr(enso_1pctS,method=pearsonr_r)
pval1 = strat_1pctS.corr(enso_1pctS,method=pearsonr_pval)

corr2=strat_4xS.corr(enso_4xS,method=pearsonr_r)
pval2 = strat_4xS.corr(enso_4xS,method=pearsonr_pval)

corr3=ssp_total_stratS.corr(diff_ssp5S,method=pearsonr_r)
pval3 = ssp_total_stratS.corr(diff_ssp5S,method=pearsonr_pval)

ax[3].set_title('1pctCO$_2$: ENSO amp. vs mean strat.',fontsize=40)
ax[3].set_xlabel('$\Delta$stratification  ($^o$C)')
ax[3].set_ylabel('$\Delta$std SST ($^o$C)')
#ax[0].set_ylim(-0.2,0.6)

ax[4].set_title('abrupt4xCO$_2$: ENSO amp. vs mean strat.',fontsize=40)
ax[4].set_xlabel('$\Delta$stratification  ($^o$C)')
ax[4].set_ylabel('$\Delta$std SST ($^o$C)')
#ax[1].set_ylim(-0.2,0.6)

ax[5].set_title('ssp585: ENSO amp. vs mean strat.',fontsize=40)
ax[5].set_xlabel('$\Delta$stratification ($^o$C)')
ax[5].set_ylabel('$\Delta$std SST ($^o$C)')
ax[5].set_ylim([-0.1,0.7])

# for i in x:
#     ax[3].scatter(diff_1pct[i],diff_1pct_pr[i],s=1000,marker=markerlist[i],c=colorlist[i])
# for i in range(0,len(x1)):
#     ax[4].scatter(diff_4x[i],diff_4x_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])
#     ax[5].scatter(diff_ssp5[i],diff_ssp5_pr[i],s=1000,marker=markerlist_subset1[i],c=colorlist_subset1[i])

# ax[3].set_title('1pctCO2 ENSO SST vs ENSO rainfall',fontsize=40)
# ax[3].set_xlabel('$\Delta$std SST')
# ax[3].set_ylabel('$\Delta$std rainfall')
# ax[3].set_ylim(0.2*10**(-5),2*10**(-5))

# ax[4].set_title('4xCO2 ENSO SST vs ENSO rainfall',fontsize=40)
# ax[4].set_xlabel('$\Delta$std SST')
# ax[4].set_ylabel('$\Delta$std rainfall')
# ax[4].set_ylim(0.2*10**(-5),2*10**(-5))

# ax[5].set_title('SSP585 ENSO SST vs ENSO rainfall',fontsize=40)
# ax[5].set_xlabel('$\Delta$std SST')
# ax[5].set_ylabel('$\Delta$std rainfall')
# ax[5].set_ylim(0.2*10**(-5),2*10**(-5))



# plt.figtext(0.05,0.92,'R='+str("%.2f" % corr1)+',  p='+str("%.2f" % pval1))
# plt.figtext(0.4,0.92,'R='+str("%.2f" % corr2)+',  p='+str("%.2f" % pval2))
# plt.figtext(0.85,0.92,'R='+str("%.2f" % corr3)+',  p='+str("%.2f" % pval3))

plt.figtext(0.05,0.40,'R='+str("%.2f" % corr1)+',  p='+str("%.2f" % pval1))
plt.figtext(0.4,0.40,'R='+str("%.2f" % corr2)+',  p='+str("%.2f" % pval2))
plt.figtext(0.75,0.40,'R='+str("%.2f" % corr3)+',  p='+str("%.2f" % pval2))

plt.figtext(0.05,0.93,'R='+str("%.2f" % corr4)+',  p='+str("%.2f" % pval4))
plt.figtext(0.4,0.93,'R='+str("%.2f" % corr5)+',  p='+str("%.2f" % pval5))
plt.figtext(0.75,0.93,'R='+str("%.2f" % corr6)+',  p='+str("%.2f" % pval6))

# plt.figtext(0.05,0.40,'R='+str("%.2f" % corr7)+',  p='+str("%.2f" % pval7))
# plt.figtext(0.4,0.40,'R='+str("%.2f" % corr8)+',  p='+str("%.2f" % pval8))
# plt.figtext(0.75,0.40,'R='+str("%.2f" % corr9)+',  p='+str("%.2f" % pval9))

