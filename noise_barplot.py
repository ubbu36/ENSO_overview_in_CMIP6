#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 22:18:06 2021

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



model_names=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']

model_names_1=['CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','HadGEM3-GC31-LL','IPSL-CM6A-LR',\
                    'MIROC6', 'MPI-ESM1-2-LR','UKESM1-0-LL']


model_names_2=['CanESM5','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','IPSL-CM6A-LR',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']



noise_C=[0]*len(model_names)
noise_4=[0]*len(model_names)
noise_1pct=[0]*len(model_names)
noise_ssp5=[0]*len(model_names_1)
noise_ssp5_max=[0]*len(model_names_1)
noise_ssp5_min=[0]*len(model_names_1)

noise_ssp1=[0]*len(model_names_2)
noise_ssp1_max=[0]*len(model_names_2)
noise_ssp1_min=[0]*len(model_names_2)

noise_4_diff=[0]*len(model_names)
noise_1pct_diff=[0]*len(model_names)
noise_ssp5_diff=[0]*len(model_names_1)
noise_C1=[0]*len(model_names_1)

for x in range(0,len(model_names)):
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names[x]+'_control.nc')
    noise=noise['__xarray_dataarray_variable__'].std('time')
    noise=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')
    noise_C[x]=noise
    noise.to_netcdf('/Users/ullaheede_1/Downloads/calc_noise_'+model_names[x]+'_control.nc')
    
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names[x]+'_4x.nc')
    noise=noise['__xarray_dataarray_variable__'].isel(time=slice(0,150*12)).std('time')
    noise=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')

    noise_4[x]=noise
    noise.to_netcdf('/Users/ullaheede_1/Downloads/calc_noise_'+model_names[x]+'_4x.nc')
    noise_4_diff[x]=(noise_4[x]-noise_C[x]).values
    
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names[x]+'_1pct.nc')
    noise=noise['__xarray_dataarray_variable__'].isel(time=slice(0,150*12)).std('time')
    noise=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')

    noise_1pct[x]=noise
    noise.to_netcdf('/Users/ullaheede_1/Downloads/calc_noise_'+model_names[x]+'_1pct.nc')
    noise_1pct_diff[x]=(noise_1pct[x]-noise_C[x]).values
    
    
    
for x in range(0,len(model_names_1)):
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names_1[x]+'_control.nc')
    noise=noise['__xarray_dataarray_variable__'].std('time')
    noise=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')
    noise_C1[x]=noise
    
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names_1[x]+'_1.nc')
    noise=noise['__xarray_dataarray_variable__'].isel(time=slice(0,85*12)).std('time')
    noise1=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')
    noise1.to_netcdf('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_1[x]+'_1.nc')
    
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names_1[x]+'_2.nc')
    noise=noise['__xarray_dataarray_variable__'].isel(time=slice(0,85*12)).std('time')
    noise2=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')
    noise2.to_netcdf('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_1[x]+'_2.nc')
    
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names_1[x]+'_3.nc')
    noise=noise['__xarray_dataarray_variable__'].isel(time=slice(0,85*12)).std('time')
    noise3=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')
    noise3.to_netcdf('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_1[x]+'_3.nc')
    
    noise_ssp5[x]=(noise1+noise2+noise3)/3
    noise_ssp5_max[x]=max(noise1,noise2,noise3)-noise_ssp5[x]
    noise_ssp5_min[x]=noise_ssp5[x]-min(noise1,noise2,noise3)
    noise_ssp5_diff[x]=(noise_ssp5[x]-noise_C1[x]).values

for x in range(0,len(model_names_2)):
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names_2[x]+'_control.nc')
    noise=noise['__xarray_dataarray_variable__'].std('time')
    noise=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')
    noise_C1[x]=noise
    
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names_2[x]+'_ssp126_1.nc')
    noise=noise['__xarray_dataarray_variable__'].isel(time=slice(0,85*12)).std('time')
    noise1=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')
    noise1.to_netcdf('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_2[x]+'_ssp126_1.nc')
    
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names_2[x]+'_ssp126_2.nc')
    noise=noise['__xarray_dataarray_variable__'].isel(time=slice(0,85*12)).std('time')
    noise2=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')
    noise2.to_netcdf('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_2[x]+'_ssp126_2.nc')
    
    noise=xr.open_dataset('/Users/ullaheede_1/Downloads/spatial_noise_'+model_names_2[x]+'_ssp126_3.nc')
    noise=noise['__xarray_dataarray_variable__'].isel(time=slice(0,85*12)).std('time')
    noise3=noise.sel(lon=slice(150,280),lat=slice(-5,5)).mean('lon').mean('lat')
    noise3.to_netcdf('/Users/ullaheede_1/Downloads/calc_noise_'+model_names_2[x]+'_ssp126_3.nc')
    
    noise_ssp1[x]=(noise1+noise2+noise3)/3
    noise_ssp1_max[x]=max(noise1,noise2,noise3)-noise_ssp1[x]
    noise_ssp1_min[x]=noise_ssp1[x]-min(noise1,noise2,noise3)

    
plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'hatch.color': '0.1'})  
x = np.arange(len(model_names))  # the label locations
x1=array([ 2,3,4,5,6,7,8,9,10,11,13,14,16,17])

    
x2=array([ 2,4,5,6,7,8,9,10,13,14,15,16,17])

width = 0.15  # the width of the bars
fig = figure(figsize=(47,25))
gs = gridspec.GridSpec(2, 3)
ax1 = plt.subplot(gs[0, 0:3])
ax2 = plt.subplot(gs[1, 0:1])
ax3 = plt.subplot(gs[1, 1:2])
ax4 = plt.subplot(gs[1, 2:3])


fig = gcf()
gs.tight_layout(fig,h_pad=12,w_pad=4.5)
ax = [ax1, ax2, ax3, ax4]

plt.figtext(0.04, 0.98, 'a)')
plt.figtext(0.04, 0.37, 'b)')
plt.figtext(0.37, 0.37, 'c)')
plt.figtext(0.70, 0.37, 'd)')


ax[0].bar(x - width*2, noise_C, width,label='PiControl',color='k')
ax[0].bar(x - width, noise_4, width,label='abrupt-4xCO2',color='#ff7f0e')
ax[0].bar(x, noise_1pct, width,label='1pctCO2',color='#2ca02c')
ax[0].bar(x1+width, noise_ssp5, width,label='ssp585',yerr=[noise_ssp5_min,noise_ssp5_max], error_kw=dict(lw=5),color='#d62728')
ax[0].bar(x2+2*width, noise_ssp1, width,label='ssp126',yerr=[noise_ssp1_min,noise_ssp1_max], error_kw=dict(lw=5),color='purple')

ax[0].bar(19 - width*2, mean(noise_C), width,color='k')
ax[0].bar(19 - width, mean(noise_4), width,color='#ff7f0e')
ax[0].bar(19, mean(noise_1pct), width,color='#2ca02c')
ax[0].bar(19+width, mean(noise_ssp5), width, error_kw=dict(lw=5),color='#d62728')
ax[0].bar(19+2*width, mean(noise_ssp1), width, error_kw=dict(lw=5),color='purple')

#ax[0].bar(x + width, ssp_total, width,yerr=ssp_max, error_kw=dict(lw=5),label='ssp585')
#ax[0].bar(x2 + width*2, ssp_total_126, width,yerr=ssp_max_126, error_kw=dict(lw=5),label='ssp126')

ax[0].set_xticks(range(len(model_names)+2))

ax[0].set_xticklabels(model_names+['']+['multi-model mean'],rotation='vertical')

ax[0].set_title('atmospheric noise in CMIP6 models',fontsize=40)

ax[0].set_ylabel('noise (N/m$^2$)')
ax[0].legend(ncol=2,fontsize=38,loc='upper left')   

#%%
    

from scipy.stats import kendalltau, pearsonr, spearmanr

def pearsonr_pval(x,y):
        return pearsonr(x,y)[1]
def pearsonr_r(x,y):
        return pearsonr(x,y)[0]


    
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

w1=120
w2=150
x=20

#grad_dif_4x=[0]*len(model_names)
#grad_dif_1pct=[0]*len(model_names)
#grad_dif_ssp5=[0]*len(model_names_1)
#grad_dif_ssp1=[0]*len(model_names_2)

enso_4x=[0]*len(model_names)
enso_control=[0]*len(model_names)
enso_ssp585_1=[0]*len(model_names_1)
enso_ssp585_2=[0]*len(model_names_1)
enso_ssp585_3=[0]*len(model_names_1)
ssp_total=[0]*len(model_names_1)
ssp_max=[0]*len(model_names)
ssp_min=[0]*len(model_names)
#ssp_total_126=[0]*len(model_names_2)
#ssp_max_126=[0]*len(model_names_2)
#ssp_min_126=[0]*len(model_names_2)
#enso_ssp126=[0]*len(model_names_2)
enso_1pct=[0]*len(model_names)

#%%
for x in range(0,len(model_names)):
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq4_'+model_names[x]+'.nc')
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

    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names[x]+'.nc')
    control=control['ts']
    control=control.assign_coords(time=list(range(len(control.time))))
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    
  #  grad_dif_4x[x]=(grad1-gradC)

    enso_ssp=(enso1)
    enso_4x[x]=enso1
    
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
    
    
    enso_control[x]=ensoC.values
    
  
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1pct_'+model_names[x]+'.nc')
    control=control['ts'].isel(time=slice(0,150*12))
    control=control.assign_coords(time=list(range(len(control.time))))
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad1=(west-east).mean('time')
    nino3=control.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=control.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso1=low_nino3.std('time')
    
    
    enso_1pct[x]=enso1
    
#    grad_dif_1pct[x]=(grad1-gradC)

for x in range(0,len(model_names_1)):    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1_'+model_names_1[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,85*12))
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
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq2_'+model_names_1[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad2=(west-east).mean('time')
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso2=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_2[x]=enso2
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq3_'+model_names_1[x]+'.nc')
    ens1=ens1['ts'].isel(time=slice(0,85*12))
    ens1=ens1.assign_coords(time=list(range(len(ens1.time))))
    east=ens1.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=ens1.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    grad3=(west-east).mean('time')
    nino3=ens1.sel(lon=slice(nino3_start,nino3_stop)).mean('lon').mean('lat')
    nino4=ens1.sel(lon=slice(nino4_start,nino4_stop)).mean('lon').mean('lat')
    low_nino3 = dsp.bandpass(nino3,1/(7*12),1/(1.5*12),dim='time')
    low_nino4 = dsp.bandpass(nino4,1/(7*12),1/(1.5*12),dim='time')
    EP=low_nino3-0.5*low_nino4
    enso3=low_nino3.std('time')


    grad_ssp=(grad1+grad2+grad3)/3
#    grad_dif_ssp5[x]=grad_ssp-gradC
    
    enso_ssp=(enso1)
    enso_ssp585_3[x]=enso3
    test=xr.concat([enso1,enso2,enso3], dim='ens')
    ssp_total[x]=(enso1+enso2+enso3)/3
    ssp_max[x]=test.max('ens')-ssp_total[x]
    ssp_min[x]=ssp_total[x]-test.min('ens')


    plt.plot(ensoC)
    plt.plot(enso_ssp)
#plt.plot(control.mean('lat').mean('time'))
#plt.plot(ens1.mean('lat').mean('time'))



markerlist=np.array(['o','v','<','>','1','2','3','4','s','p','P','*','h','+','X','D','|','_','o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_'])
colorlist=np.array(['grey','brown','olive','green','cyan','blue','purple','pink','red','k','grey','brown','orange','olive','cyan','blue','purple','pink','red','k',\
    'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k', 'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k'])

 
x1=array([2, 3,  4,  5,  6,  7,  8,  9, 10, 11, 13, 14, 16, 17]).astype(int)

markerlist_subset1=markerlist[x1]
colorlist_subset1=colorlist[x1]

enso_control=np.array(enso_control)
enso_control_1=enso_control[x1]

diff_1pct=enso_1pct-enso_control
diff_ssp5=ssp_total-enso_control_1
diff_4x = enso_4x-enso_control

diff_1pctS=pd.Series(diff_1pct)
diff_ssp5S=pd.Series(diff_ssp5)
diff_4xS=pd.Series(diff_4x)


noise_4_diffS=pd.Series(noise_4_diff)
noise_1pct_diffS=pd.Series(noise_1pct_diff)
noise_ssp5_diffS=pd.Series(noise_ssp5_diff)

corr1=diff_1pctS.corr(noise_1pct_diffS,method=pearsonr_r)
pval1 = diff_1pctS.corr(noise_1pct_diffS,method=pearsonr_pval)

corr2=diff_4xS.corr(noise_4_diffS,method=pearsonr_r)
pval2 = diff_4xS.corr(noise_4_diffS,method=pearsonr_pval)

corr3=diff_ssp5S.corr(noise_ssp5_diffS,method=pearsonr_r)
pval3 = diff_ssp5S.corr(noise_ssp5_diffS,method=pearsonr_pval)

x = np.arange(len(model_names)) 
x1 = np.arange(len(model_names_1)) 
for i in x:
    ax[1].scatter(diff_1pctS[i],noise_1pct_diffS[i],s=800,marker=markerlist[i],c=colorlist[i])
    ax[2].scatter(diff_4xS[i],noise_4_diffS[i],s=800,marker=markerlist[i],c=colorlist[i])
for i in x1:
    ax[3].scatter(diff_ssp5S[i],noise_ssp5_diffS[i],s=800,marker=markerlist_subset1[i],c=colorlist_subset1[i])


ax[1].set_title('1pctCO2 ENSO SST vs noise',fontsize=40)
ax[1].set_xlabel('$\Delta$ENSO SST amp. ($^o$C)')
ax[1].set_ylabel('$\Delta$noise (N/m$^2$)')
ax[1].set_ylim(-2.5*10**(-3),4.5*10**(-3))
ax[1].set_xlim(-0.2,0.5)

ax[2].set_title('abrupt4xCO2 ENSO SST vs noise',fontsize=40)
ax[2].set_xlabel('$\Delta$ENSO SST amp. ($^o$C)')
ax[2].set_ylabel('$\Delta$noise (N/m$^2$)')
ax[2].set_ylim(-2.5*10**(-3),4.5*10**(-3))
ax[2].set_xlim(-0.2,0.5)

ax[3].set_title('SSP585 ENSO SST vs noise',fontsize=40)
ax[3].set_xlabel('$\Delta$ENSO SST amp. ($^o$C)')
ax[3].set_ylabel('$\Delta$noise (N/m$^2$)')
ax[3].set_ylim(-2.5*10**(-3),4.5*10**(-3))
ax[3].set_xlim(-0.2,0.5)

plt.figtext(0.05,0.31,'R='+str("%.2f" % corr1)+',  p='+str("%.2f" % pval1))
plt.figtext(0.39,0.31,'R='+str("%.2f" % corr2)+',  p='+str("%.2f" % pval2))
plt.figtext(0.74,0.31,'R='+str("%.2f" % corr3)+',  p='+str("%.2f" % pval3))


#%%




