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
from scipy.stats import kendalltau, pearsonr, spearmanr

def pearsonr_pval(x,y):
        return pearsonr(x,y)[1]
def pearsonr_r(x,y):
        return pearsonr(x,y)[0]

model_names=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5-CanOE','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']

model_names_1=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2-WACCM','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IPSL-CM6A-LR','KACE-1-0-G',\
                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']


#model_names_1=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2-WACCM','CESM2',\
#             'EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
#                 'FIO-ESM-2-0','HadGEM3-GC31-LL','HadGEM3-GC31-MM','KACE-1-0-G',\
#                    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-LR','UKESM1-0-LL']    

model_names_2=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5-CanOE','CanESM5','CESM2',\
             'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-g3',\
                 'FIO-ESM-2-0','IPSL-CM6A-LR','KACE-1-0-G',\
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

w1=120
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

    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_1[x]+'.nc')
    control=control['ts']
    control=control.assign_coords(time=list(range(len(control.time))))
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    
    grad_dif_4x[x]=(grad1.values-gradC.values)

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
    
    
    enso_control[x]=ensoC
    
  
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1pct_'+model_names[x]+'.nc')
    control=control['ts'].isel(time=slice(50*12,150*12))
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
    
    grad_dif_1pct[x]=(grad1.values-gradC.values)
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1_'+model_names[x]+'.nc')
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


    grad_ssp=grad1.values

    enso_ssp=(enso1)
    enso_ssp585_1[x]=enso1
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq2_'+model_names[x]+'.nc')
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


    grad_ssp=grad1.values

    enso_ssp=(enso1)
    enso_ssp585_2[x]=enso2
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq3_'+model_names[x]+'.nc')
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
    grad_dif_ssp5[x]=(grad_ssp.values-gradC.values)
    
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

#%%
for x in range(0,len(model_names_2)):
    
    control=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eqC_'+model_names_2[x]+'.nc')
    control=control['ts']
    control=control.assign_coords(time=list(range(len(control.time))))
    east=control.sel(lon=slice(e1,e2)).mean('lon').mean('lat')
    west=control.sel(lon=slice(w1,w2)).mean('lon').mean('lat')
    gradC=(west-east).mean('time')
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq1ssp126_'+model_names_2[x]+'.nc')
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
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq2ssp126_'+model_names_2[x]+'.nc')
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
    enso2=low_nino3.std('time')


    grad_ssp=(grad1)

    enso_ssp=(enso1)
    enso_ssp585_2[x]=enso2
    
    ens1=xr.open_dataset('/Users/ullaheede_1/Downloads/ts_eq3ssp126_'+model_names_2[x]+'.nc')
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


    grad_ssp=(grad1.values+grad2.values+grad3.values)/3
    grad_dif_ssp1[x]=grad_ssp-gradC
    enso_ssp=(enso1)
    enso_ssp585_3[x]=enso3
    test=xr.concat([enso1,enso2,enso3], dim='ens')
    ssp_total_126[x]=test.mean('ens')
    ssp_max_126[x]=test.max('ens')-ssp_total_126[x]
    ssp_min_126[x]=ssp_total_126[x]-test.min('ens')

#%%

markerlist=np.array(['o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_','o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_'])
colorlist=np.array(['grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k','grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k',\
    'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k', 'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k'])

plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'hatch.color': '0.1'})  
x = np.arange(len(model_names)).astype(int)  # the label locations

x1=array([ 0,  1, 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,16,17,18,19])

#x1=array([ 0,  1, 3,  4,    7,  8,  9, 11, 12, 13, 14, 15,16,17,18,19])
x1 = x1.astype(int)

x2=array([ 0,  1, 2, 3,  5,  6,  7,  8,  9, 10,11, 14, 15,16,17,18,19])

width = 0.15  # the width of the bars
fig = figure(figsize=(47,25))
gs = gridspec.GridSpec(2, 3)
ax1 = plt.subplot(gs[0, 0:3])
ax2 = plt.subplot(gs[1, 0:1])

ax3 = plt.subplot(gs[1, 1:2])
ax4 = plt.subplot(gs[1, 2:3])

fig = gcf()
gs.tight_layout(fig,h_pad=12,w_pad=1.5)
ax = [ax1, ax2, ax3, ax4]

plt.figtext(0.04, 0.98, 'a)')
plt.figtext(0.04, 0.37, 'b)')
plt.figtext(0.37, 0.37, 'c)')
plt.figtext(0.70, 0.37, 'd)')


ax[0].bar(x - width*2, enso_control, width,label='PiControl',color='k')
ax[0].bar(x1 - width, enso_4x, width,label='abrupt-4xCO2',color='orange')
ax[0].bar(x, enso_1pct, width,label='1pctCO2',color='green')
ax[0].bar(x + width, ssp_total, width,yerr=[ssp_min,ssp_max], error_kw=dict(lw=5),label='ssp585',color='red')
ax[0].bar(x2 + width*2, ssp_total_126, width,yerr=[ssp_min_126,ssp_max_126], error_kw=dict(lw=5),label='ssp126',color='purple')

ax[0].bar(21 - width*2,mean(enso_control),width, color='k')
ax[0].bar(21 - width,mean(enso_4x),width, color='orange')
ax[0].bar(21,mean(enso_1pct),width, color='green')
ax[0].bar(21 + width,mean(ssp_total),width, color='red')
ax[0].bar(21 + width*2,mean(ssp_total_126),width, color='purple')

model_names3=model_names+['']+['multi-model mean']

ax[0].set_xticks(range(len(model_names3)))

ax[0].set_xticklabels(model_names3,rotation='vertical')

ax[0].set_title('ENSO in CMIP6 models',fontsize=40)

ax[0].set_ylabel('ENSO SST amp. ($^o$C)')
ax[0].legend(ncol=2,fontsize=38)

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

diff_ssp5_subset1=diff_ssp5[x1]
diff_ssp5_subset2=diff_ssp5[x2]
markerlist_subset1=markerlist[x1]
markerlist_subset2=markerlist[x2]
colorlist_subset1=colorlist[x1]
colorlist_subset2=colorlist[x2]
for i in x:
    ax[1].scatter(diff_ssp5[i],diff_1pct[i], marker=markerlist[i],s=800,c=colorlist[i])
for i in range(0,len(x1)):
    ax[2].scatter(diff_ssp5_subset1[i],diff_4x[i],marker=markerlist_subset1[i],s=800,c=colorlist_subset1[i])
for i in range(0,len(x2)): 
    ax[3].scatter(diff_ssp5_subset2[i],diff_ssp1[i], marker=markerlist_subset2[i],s=800,c=colorlist_subset2[i])


ax[1].set_title('1pctCO2 vs ssp585',fontsize=40)
ax[1].set_ylabel('$\Delta$ENSO 1pctCO2')
ax[1].set_xlabel('$\Delta$ENSO ssp585')

ax[2].set_title('4xCO2 vs ssp585',fontsize=40)
ax[2].set_ylabel('$\Delta$ENSO 4xCO2')
ax[2].set_xlabel('$\Delta$ENSO ssp585')

ax[3].set_title('ssp585 vs ssp126',fontsize=40)
ax[3].set_xlabel('$\Delta$ENSO ssp585')
ax[3].set_ylabel('$\Delta$ENSO ssp126')
#ax[3].set_ylim(-0,1,0.7)
#ax[3].set_xlim(-0,1,0.7)

diff_1pctS=pd.Series(diff_1pct)
diff_ssp5S=pd.Series(diff_ssp5)
diff_ssp5S_subset1=pd.Series(diff_ssp5[x1])
diff_ssp5S_subset2=pd.Series(diff_ssp5[x2])

diff_4xS=pd.Series(diff_4x)
diff_ssp1S=pd.Series(diff_ssp1)

corr1=diff_1pctS.corr(diff_ssp5S,method=pearsonr_r)
pval1 = diff_1pctS.corr(diff_ssp5S,method=pearsonr_pval)

corr2=diff_4xS.corr(diff_ssp5S_subset1,method=pearsonr_r)
pval2=diff_4xS.corr(diff_ssp5S_subset1,method=pearsonr_pval)

corr3=diff_ssp1S.corr(diff_ssp5S_subset2,method=pearsonr_r)
pval3=diff_ssp1S.corr(diff_ssp5S_subset2,method=pearsonr_pval)

plt.figtext(0.05,0.31,'R='+str("%.2f" % corr1)+',  p='+str("%.2f" % pval1))
plt.figtext(0.38,0.31,'R='+str("%.2f" % corr2)+',  p='+str("%.2f" % pval2))
plt.figtext(0.71,0.31,'R='+str("%.2f" % corr3)+',  p='+str("%.2f" % pval3))

#%%

#markerlist=np.array(['o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_','o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_'])
#colorlist=np.array(['grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k','grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k',\
#    'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k', 'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k'])

#plt.scatter(1,1,marker=markerlist[19],s=2000,c=colorlist[19])

