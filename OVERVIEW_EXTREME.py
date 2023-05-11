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

model_names_1=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']
    
model_names_2=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5-CanOE','CanESM5','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','IPSL-CM6A-LR','KACE-1-0-G',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']
    
#model_names=['CanESM5','CESM2-WACCM','CESM2',\
#             'CNRM-CM6-1','CNRM-ESM2-1','FGOALS-g3',\
#                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
#                     'MPI-ESM1-2-LR','UKESM1-0-LL']


nino3_start=210
nino3_stop=280

nino4_start=160
nino4_stop=210
    
e1=200
e2=280

w1=120
w2=180
x=20

grad_dif_4x=[0]*len(model_names)
grad_dif_4x=[0]*len(model_names)
enso_4x_ex=[0]*len(model_names_1)
enso_control_ex=[0]*len(model_names)
enso_ssp585_1=[0]*len(model_names)
enso_ssp585_2=[0]*len(model_names)
enso_ssp585_3=[0]*len(model_names)
ssp_total_ex=[0]*len(model_names)
ssp_max_ex=[0]*len(model_names)
ssp_min_ex=[0]*len(model_names)
ssp_total_126_ex=[0]*len(model_names_2)
ssp_max_126_ex=[0]*len(model_names_2)
ssp_min_126_ex=[0]*len(model_names_2)
enso_ssp126_1=[0]*len(model_names_2)
enso_ssp126_2=[0]*len(model_names_2)
enso_ssp126_3=[0]*len(model_names_2)
enso_1pct_ex=[0]*len(model_names)

#%%
for x in range(0,len(model_names_1)):
    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_1[x]+'.nc')
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
    c=ensoC.values
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq4_'+model_names_1[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,150*12))
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

    peaks_4=scipy.signal.find_peaks(low_nino3, height=c*2, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    

    enso_4x_ex[x]=len(peaks_4[0])/len(ens1)*10*12
    
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
    c=ensoC.values
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,150*12))
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

    peaks_C=scipy.signal.find_peaks(low_nino3, height=c*2, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    

    enso_control_ex[x]=len(peaks_C[0])/len(ens1)*10*12
    
 ####################################################   
    
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
    c=ensoC.values
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1pct_'+model_names[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,150*12))
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

    peaks_1pct=scipy.signal.find_peaks(low_nino3, height=c*2, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    

    enso_1pct_ex[x]=len(peaks_1pct[0])/len(ens1)*10*12
####################################################    

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
    c=ensoC.values
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1_'+model_names[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,150*12))
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

    peaks_1=scipy.signal.find_peaks(low_nino3, height=c*2, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    

    enso_ssp585_1[x]=len(peaks_1[0])/len(ens1)*10*12
    
 ####################################################     
    
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
    c=ensoC.values
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq2_'+model_names[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,150*12))
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

    peaks_2=scipy.signal.find_peaks(low_nino3, height=c*2, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    

    enso_ssp585_2[x]=len(peaks_2[0])/len(ens1)*10*12
    
####################################################     
    
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
    c=ensoC.values
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq3_'+model_names[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,150*12))
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

    peaks_3=scipy.signal.find_peaks(low_nino3, height=c*2, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    

    enso_ssp585_3[x]=len(peaks_3[0])/len(ens1)*10*12

    ssp_total_ex[x]=(enso_ssp585_1[x]+enso_ssp585_2[x]+enso_ssp585_3[x])/3
    ssp_max_ex[x]=max([enso_ssp585_1[x],enso_ssp585_2[x],enso_ssp585_3[x]])-ssp_total_ex[x]
    ssp_min_ex[x]=ssp_total_ex[x]-min([enso_ssp585_1[x],enso_ssp585_2[x],enso_ssp585_3[x]])
    
for x in range(0,len(model_names_2)):
####################################################    

    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_2[x]+'.nc')
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
    c=ensoC.values
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1ssp126_'+model_names_2[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,150*12))
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

    peaks_1=scipy.signal.find_peaks(low_nino3, height=c*2, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    

    enso_ssp126_1[x]=len(peaks_1[0])/len(ens1)*10*12
    
 ####################################################     
    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_2[x]+'.nc')
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
    c=ensoC.values
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq2ssp126_'+model_names_2[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,150*12))
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

    peaks_2=scipy.signal.find_peaks(low_nino3, height=c*2, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    

    enso_ssp126_2[x]=len(peaks_2[0])/len(ens1)*10*12
    
####################################################     
    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_2[x]+'.nc')
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
    c=ensoC.values
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq3ssp126_'+model_names_2[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,150*12))
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

    peaks_3=scipy.signal.find_peaks(low_nino3, height=c*2, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    

    enso_ssp126_3[x]=len(peaks_3[0])/len(ens1)*10*12
    
    ssp_total_126_ex[x]=(enso_ssp126_1[x]+enso_ssp126_2[x]+enso_ssp126_3[x])/3
    ssp_max_126_ex[x]=max([enso_ssp126_1[x],enso_ssp126_2[x],enso_ssp126_3[x]])-ssp_total_126_ex[x]
    ssp_min_126_ex[x]=ssp_total_126_ex[x]-min([enso_ssp126_1[x],enso_ssp126_2[x],enso_ssp126_3[x]])
#%%

markerlist=np.array(['o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_','o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_'])
colorlist=np.array(['grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k','grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k',\
    'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k', 'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k'])

plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'hatch.color': '0.1'})  
x = np.arange(len(model_names))  # the label locations
x1=array([ 0,  1, 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,16,17,18,19])
x2=array([ 0,  1, 2, 3,  5,  6,  7,  8,  9, 10,11, 14, 15,16,17,18,19])
width = 0.15  # the width of the bars
fig = figure(figsize=(47,25))
gs = gridspec.GridSpec(2, 3)
ax1 = plt.subplot(gs[0, 0:3])
ax2 = plt.subplot(gs[1, 0:1])
ax3 = plt.subplot(gs[1, 1:2])
ax4 = plt.subplot(gs[1, 2:3])

markerlist_subset1=markerlist[x1]
markerlist_subset2=markerlist[x2]
colorlist_subset1=colorlist[x1]
colorlist_subset2=colorlist[x2]

fig = gcf()
gs.tight_layout(fig,h_pad=12,w_pad=2)
ax = [ax1, ax2, ax3, ax4]

plt.figtext(0.04, 0.98, 'a)')
plt.figtext(0.04, 0.37, 'b)')
plt.figtext(0.37, 0.37, 'c)')
plt.figtext(0.70, 0.37, 'd)')

ax[0].bar(x - width*2, enso_control_ex, width,label='PiControl',color='k')
ax[0].bar(x1 - width, enso_4x_ex, width,label='abrupt-4xCO2',color='orange')
ax[0].bar(x, enso_1pct_ex, width,label='1pctCO2',color='green')
ax[0].bar(x + width, ssp_total_ex, width,yerr=[ssp_min_ex,ssp_max_ex], error_kw=dict(lw=5),label='ssp585',color='red')
ax[0].bar(x2 + width*2, ssp_total_126_ex, width,yerr=[ssp_min_126_ex,ssp_max_126_ex], error_kw=dict(lw=5),label='ssp126',color='purple')

ax[0].bar(21 - width*2,mean(enso_control_ex),width, color='k')
ax[0].bar(21 - width,mean(enso_4x_ex),width, color='orange')
ax[0].bar(21,mean(enso_1pct_ex),width, color='green')
ax[0].bar(21 + width,mean(ssp_total_ex),width, color='red')
ax[0].bar(21 + width*2,mean(ssp_total_126_ex),width, color='purple')

model_names3=model_names+['']+['multi-model mean']

ax[0].set_xticks(range(len(model_names3)))

ax[0].set_xticklabels(model_names3,rotation='vertical')

ax[0].set_title('extreme El Niño events',fontsize=40)

ax[0].set_ylabel('ex events per dec')
ax[0].legend(ncol=2,fontsize=38)



enso_control_ex=np.array(enso_control_ex)

enso_control_1_ex=enso_control_ex[x1]
enso_control_2_ex=enso_control_ex[x2]

enso_4x_ex=np.array(enso_4x_ex)
enso_1pct_ex=np.array(enso_1pct_ex)


ssp_total_ex=np.array(ssp_total_ex)
ssp_total_126_ex=np.array(ssp_total_126_ex)

diff_4x_ex = (enso_4x_ex-enso_control_1_ex)/enso_control_1_ex*100
 
diff_1pct_ex=(enso_1pct_ex-enso_control_ex)/enso_control_ex*100
diff_ssp5_ex=(ssp_total_ex-enso_control_ex)/enso_control_ex*100
diff_ssp1_ex=(ssp_total_126_ex-enso_control_2_ex)/enso_control_2_ex*100

diff_1pctS=pd.Series(diff_1pct_p)
diff_ssp5S=pd.Series(diff_ssp5_p)
diff_4xS=pd.Series(diff_4x_p)

diff_1pctS_ex=pd.Series(diff_1pct_ex)
diff_ssp5S_ex=pd.Series(diff_ssp5_ex)
diff_4xS_ex=pd.Series(diff_4x_ex)

corr1=diff_1pctS.corr(diff_1pctS_ex,method=pearsonr_r)
pval1 = diff_1pctS.corr(diff_1pctS_ex,method=pearsonr_pval)

corr2=diff_4xS.corr(diff_4xS_ex,method=pearsonr_r)
pval2=diff_4xS.corr(diff_4xS_ex,method=pearsonr_pval)

corr3=diff_ssp5S.corr(diff_ssp5S_ex,method=pearsonr_r)
pval3=diff_ssp5S.corr(diff_ssp5S_ex,method=pearsonr_pval)

plt.figtext(0.05,0.31,'R='+str("%.2f" % corr1)+',  p='+str("%.2f" % pval1))
plt.figtext(0.38,0.31,'R='+str("%.2f" % corr2)+',  p='+str("%.2f" % pval2))
plt.figtext(0.71,0.31,'R='+str("%.2f" % corr3)+',  p='+str("%.2f" % pval3))

#ax[1].scatter(diff_1pct_p,diff_1pct_ex,s=400)
#ax[2].scatter(diff_4x_p,diff_4x_ex,s=400)
#ax[3].scatter(diff_ssp5_p,diff_ssp5_ex,s=400)

for i in x:
    ax[1].scatter(diff_1pct_p[i],diff_1pct_ex[i], marker=markerlist[i],s=800,c=colorlist[i])
for i in range(0,len(x1)):
    ax[2].scatter(diff_4x_p[i],diff_4x_ex[i],marker=markerlist_subset1[i],s=800,c=colorlist_subset1[i])
for i in range(0,len(x)): 
    ax[3].scatter(diff_ssp5_p[i],diff_ssp5_ex[i], marker=markerlist[i],s=800,c=colorlist[i])

ax[1].set_title('1pctCO2',fontsize=40)
ax[1].set_xlabel('% change in ENSO amp.')
ax[1].set_ylabel('% change in ex. events')
ax[1].set_ylim(-150,320)
ax[1].set_xlim(-25,90)

ax[2].set_title('4xCO2',fontsize=40)
ax[2].set_xlabel('% change in in ENSO amp.')
ax[2].set_ylabel('% change in ex. events')
ax[2].set_ylim(-150,320)
ax[2].set_xlim(-25,90)

ax[3].set_title('ssp585',fontsize=40)
ax[3].set_xlabel('% change in ENSO amp.')
ax[3].set_ylabel('% change in ex. events')
ax[3].set_ylim(-150,320)
ax[3].set_xlim(-25,90)

#%%
data = pd.DataFrame({'x': grad_dif_4x, 'y': diff_4x_p, 'z': diff_4x_ex})
#data = pd.DataFrame({'y': diff_4x_p, 'z': diff_4x_ex})
from statsmodels.formula.api import ols
#model = ols("z ~ y", data).fit()

model = ols("z ~ y + x", data).fit()

# Print the summary
print(model.summary())