#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 13:57:50 2022

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
import xrscipy.signal as dsp
from pylab import *
import matplotlib.gridspec as gridspec

plt.rcParams.update({'font.size': 22})

fig = figure(figsize=(24,35))
gs = gridspec.GridSpec(7, 3)
ax1 = plt.subplot(gs[0, 0:1])
ax2 = plt.subplot(gs[0, 1:2])
ax3 = plt.subplot(gs[0, 2:3])

ax4 = plt.subplot(gs[1, 0:1])
ax5 = plt.subplot(gs[1, 1:2])
ax6 = plt.subplot(gs[1, 2:3])

ax7 = plt.subplot(gs[2, 0:1])
ax8 = plt.subplot(gs[2, 1:2])
ax9 = plt.subplot(gs[2, 2:3])

ax10 = plt.subplot(gs[3, 0:1])
ax11 = plt.subplot(gs[3, 1:2])
ax12 = plt.subplot(gs[3, 2:3])

ax13 = plt.subplot(gs[4, 0:1])
ax14 = plt.subplot(gs[4, 1:2])
ax15 = plt.subplot(gs[4, 2:3])

ax16 = plt.subplot(gs[5, 0:1])
ax17 = plt.subplot(gs[5, 1:2])
ax18 = plt.subplot(gs[5, 2:3])

ax19 = plt.subplot(gs[6, 0:1])
ax20 = plt.subplot(gs[6, 1:2])
ax21 = plt.subplot(gs[6, 2:3])


#ax6 = plt.subplot(gs[0, 5:6])

fig = gcf()
gs.tight_layout(fig,h_pad=2.2,w_pad=2)
ax = [ax1,ax2,ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17, ax18, ax19, ax20, ax21]

plt.figtext(0.01, 0.99, 'CanESM5',fontweight='bold')
plt.figtext(0.01, 0.846, 'CESM2',fontweight='bold')
plt.figtext(0.01, 0.703, 'CNRM-CM6-1',fontweight='bold')
plt.figtext(0.01, 0.56, 'HadGEM3-GC31-LL',fontweight='bold')
plt.figtext(0.01, 0.416, 'IPSL-CM6A-LR',fontweight='bold')
plt.figtext(0.01, 0.275, 'MIROC6',fontweight='bold')
plt.figtext(0.01, 0.132, 'MIROC-ES2L',fontweight='bold')


ax[7].set_xticks([])
ax[7].set_yticks([])

ax[13].set_xticks([])
ax[13].set_yticks([])


yr='4x.nc'
yr1='4x'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/nino_temp_control.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/eq_zo_wind_control.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[0].scatter(nino_temp1,nino_Q1)
ax[0].scatter(nino_temp2,nino_Q2)
ax[0].scatter(nino_temp3,nino_Q3)
ax[0].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[0].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[0].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[0].set_xlabel('SST anom. ($^o$C)')
ax[0].set_ylabel('N/m$^2$')
ax[0].set_title('control')
ax[0].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)

ax[0].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[0].yaxis.set_label_coords(-.15, .5)
nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)

from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[1].scatter(nino_temp1,nino_Q1)
ax[1].scatter(nino_temp2,nino_Q2)
ax[1].scatter(nino_temp3,nino_Q3)
ax[1].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[1].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[1].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[1].set_xlabel('SST anom. ($^o$C)')
ax[1].set_ylabel('N/m$^2$')
ax[1].set_title('4xCO$_2$')
ax[1].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)

ax[1].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[1].yaxis.set_label_coords(-.15, .5)

yr='ens1.nc'
yr1='ens1'
yr2='ens1'
yr3='ens1'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens1= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/nino_temp_'+yr2+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens2= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/nino_temp_'+yr3+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens3= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')


nino_temp4=xr.concat([nino_temp_ens1,nino_temp_ens2,nino_temp_ens3],'time')

nino_temp1=nino_temp4.where(nino_temp4>1)
nino_temp2=nino_temp4.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_temp3=nino_temp4.where(nino_temp4<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens1= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/eq_zo_wind_'+yr2+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens2= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CanESM5/eq_zo_wind_'+yr3+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens3= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.concat([nino_Q_ens1,nino_Q_ens2,nino_Q_ens3],'time')


nino_Q1=nino_Q.where(nino_temp4>1)
nino_Q2=nino_Q.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_Q3=nino_Q.where(nino_temp4<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[2].scatter(nino_temp1,nino_Q1)
ax[2].scatter(nino_temp2,nino_Q2)
ax[2].scatter(nino_temp3,nino_Q3)
ax[2].plot(nino_temp4,(nino_temp1*slope1+intercept1),'k')
ax[2].plot(nino_temp4,(nino_temp2*slope2+intercept2),'k')
ax[2].plot(nino_temp4,(nino_temp3*slope3+intercept3),'k')
ax[2].set_xlabel('SST anom. ($^o$C)')
ax[2].set_ylabel('N/m$^2$')
ax[2].set_title('SSP585 (3 ens. mem. )')
ax[2].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
ax[2].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
ax[2].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

ax[2].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
ax[2].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
ax[2].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
ax[2].yaxis.set_label_coords(-.15, .5)
#%%
yr='4x.nc'
yr1='4x'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/nino_temp_control.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/eq_zo_wind_control.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[3].scatter(nino_temp1,nino_Q1)
ax[3].scatter(nino_temp2,nino_Q2)
ax[3].scatter(nino_temp3,nino_Q3)
ax[3].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[3].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[3].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[3].set_xlabel('SST anom. ($^o$C)')
ax[3].set_ylabel('N/m$^2$')
ax[3].set_title('control')
ax[3].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[3].transAxes)
ax[3].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[3].transAxes)
ax[3].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[3].transAxes)

ax[3].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[3].transAxes)
ax[3].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[3].transAxes)
ax[3].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[3].transAxes)
ax[3].yaxis.set_label_coords(-.15, .5)

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)

from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[4].scatter(nino_temp1,nino_Q1)
ax[4].scatter(nino_temp2,nino_Q2)
ax[4].scatter(nino_temp3,nino_Q3)
ax[4].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[4].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[4].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[4].set_xlabel('SST anom. ($^o$C)')
ax[4].set_ylabel('N/m$^2$')
ax[4].set_title('4xCO$_2$')
ax[4].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[4].transAxes)
ax[4].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[4].transAxes)
ax[4].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[4].transAxes)

ax[4].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[4].transAxes)
ax[4].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[4].transAxes)
ax[4].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[4].transAxes)
ax[4].yaxis.set_label_coords(-.15, .5)

yr='ens1.nc'
yr1='ens1'
yr2='ens1'
yr3='ens1'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens1= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/nino_temp_'+yr2+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens2= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/nino_temp_'+yr3+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens3= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')


nino_temp4=xr.concat([nino_temp_ens1,nino_temp_ens2,nino_temp_ens3],'time')

nino_temp1=nino_temp4.where(nino_temp4>1)
nino_temp2=nino_temp4.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_temp3=nino_temp4.where(nino_temp4<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens1= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/eq_zo_wind_'+yr2+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens2= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CESM2/eq_zo_wind_'+yr3+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens3= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.concat([nino_Q_ens1,nino_Q_ens2,nino_Q_ens3],'time')


nino_Q1=nino_Q.where(nino_temp4>1)
nino_Q2=nino_Q.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_Q3=nino_Q.where(nino_temp4<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[5].scatter(nino_temp1,nino_Q1)
ax[5].scatter(nino_temp2,nino_Q2)
ax[5].scatter(nino_temp3,nino_Q3)
ax[5].plot(nino_temp4,(nino_temp1*slope1+intercept1),'k')
ax[5].plot(nino_temp4,(nino_temp2*slope2+intercept2),'k')
ax[5].plot(nino_temp4,(nino_temp3*slope3+intercept3),'k')
ax[5].set_xlabel('SST anom. ($^o$C)')
ax[5].set_ylabel('N/m$^2$')
ax[5].set_title('SSP585 (3 ens. mem. )')
ax[5].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[5].transAxes)
ax[5].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[5].transAxes)
ax[5].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[5].transAxes)

ax[5].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[5].transAxes)
ax[5].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[5].transAxes)
ax[5].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[5].transAxes)
ax[5].yaxis.set_label_coords(-.15, .5)
#%%
yr='4x.nc'
yr1='4x'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CNRM-CM6-1/nino_temp_control.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CNRM-CM6-1/eq_zo_wind_control.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[6].scatter(nino_temp1,nino_Q1)
ax[6].scatter(nino_temp2,nino_Q2)
ax[6].scatter(nino_temp3,nino_Q3)
ax[6].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[6].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[6].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[6].set_xlabel('SST anom. ($^o$C)')
ax[6].set_ylabel('N/m$^2$')
ax[6].set_title('control')
ax[6].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[6].transAxes)
ax[6].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[6].transAxes)
ax[6].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[6].transAxes)

ax[6].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[6].transAxes)
ax[6].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[6].transAxes)
ax[6].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[6].transAxes)
ax[6].yaxis.set_label_coords(-.15, .5)

yr='ens1.nc'
yr1='ens1'
yr2='ens1'
yr3='ens1'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CNRM-CM6-1/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens1= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CNRM-CM6-1/nino_temp_'+yr2+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens2= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CNRM-CM6-1/nino_temp_'+yr3+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens3= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')


nino_temp4=xr.concat([nino_temp_ens1,nino_temp_ens2,nino_temp_ens3],'time')

nino_temp1=nino_temp4.where(nino_temp4>1)
nino_temp2=nino_temp4.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_temp3=nino_temp4.where(nino_temp4<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CNRM-CM6-1/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens1= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CNRM-CM6-1/eq_zo_wind_'+yr2+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens2= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/CNRM-CM6-1/eq_zo_wind_'+yr3+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens3= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.concat([nino_Q_ens1,nino_Q_ens2,nino_Q_ens3],'time')


nino_Q1=nino_Q.where(nino_temp4>1)
nino_Q2=nino_Q.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_Q3=nino_Q.where(nino_temp4<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[8].scatter(nino_temp1,nino_Q1)
ax[8].scatter(nino_temp2,nino_Q2)
ax[8].scatter(nino_temp3,nino_Q3)
ax[8].plot(nino_temp4,(nino_temp1*slope1+intercept1),'k')
ax[8].plot(nino_temp4,(nino_temp2*slope2+intercept2),'k')
ax[8].plot(nino_temp4,(nino_temp3*slope3+intercept3),'k')
ax[8].set_xlabel('SST anom. ($^o$C)')
ax[8].set_ylabel('N/m$^2$')
ax[8].set_title('SSP585 (3 ens. mem. )')
ax[8].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[8].transAxes)
ax[8].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[8].transAxes)
ax[8].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[8].transAxes)

ax[8].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[8].transAxes)
ax[8].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[8].transAxes)
ax[8].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[8].transAxes)
ax[8].yaxis.set_label_coords(-.15, .5)
yr='4x.nc'
yr1='4x'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/nino_temp_control.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/eq_zo_wind_control.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[9].scatter(nino_temp1,nino_Q1)
ax[9].scatter(nino_temp2,nino_Q2)
ax[9].scatter(nino_temp3,nino_Q3)
ax[9].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[9].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[9].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[9].set_xlabel('SST anom. ($^o$C)')
ax[9].set_ylabel('N/m$^2$')
ax[9].set_title('control')
ax[9].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[9].transAxes)
ax[9].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[9].transAxes)
ax[9].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[9].transAxes)

ax[9].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[9].transAxes)
ax[9].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[9].transAxes)
ax[9].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[9].transAxes)
ax[9].yaxis.set_label_coords(-.15, .5)

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)

from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[10].scatter(nino_temp1,nino_Q1)
ax[10].scatter(nino_temp2,nino_Q2)
ax[10].scatter(nino_temp3,nino_Q3)
ax[10].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[10].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[10].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[10].set_xlabel('SST anom. ($^o$C)')
ax[10].set_ylabel('N/m$^2$')
ax[10].set_title('4xCO$_2$')
ax[10].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[10].transAxes)
ax[10].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[10].transAxes)
ax[10].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[10].transAxes)

ax[10].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[10].transAxes)
ax[10].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[10].transAxes)
ax[10].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[10].transAxes)
ax[10].yaxis.set_label_coords(-.15, .5)

yr='ens1.nc'
yr1='ens1'
yr2='ens1'
yr3='ens1'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens1= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/nino_temp_'+yr2+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens2= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/nino_temp_'+yr3+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens3= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')


nino_temp4=xr.concat([nino_temp_ens1,nino_temp_ens2,nino_temp_ens3],'time')

nino_temp1=nino_temp4.where(nino_temp4>1)
nino_temp2=nino_temp4.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_temp3=nino_temp4.where(nino_temp4<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens1= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/eq_zo_wind_'+yr2+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens2= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/HadGEM3-GC31-LL/eq_zo_wind_'+yr3+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens3= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.concat([nino_Q_ens1,nino_Q_ens2,nino_Q_ens3],'time')


nino_Q1=nino_Q.where(nino_temp4>1)
nino_Q2=nino_Q.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_Q3=nino_Q.where(nino_temp4<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[11].scatter(nino_temp1,nino_Q1)
ax[11].scatter(nino_temp2,nino_Q2)
ax[11].scatter(nino_temp3,nino_Q3)
ax[11].plot(nino_temp4,(nino_temp1*slope1+intercept1),'k')
ax[11].plot(nino_temp4,(nino_temp2*slope2+intercept2),'k')
ax[11].plot(nino_temp4,(nino_temp3*slope3+intercept3),'k')
ax[11].set_xlabel('SST anom. ($^o$C)')
ax[11].set_ylabel('N/m$^2$')
ax[11].set_title('SSP585 (3 ens. mem. )')
ax[11].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[11].transAxes)
ax[11].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[11].transAxes)
ax[11].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[11].transAxes)

ax[11].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[11].transAxes)
ax[11].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[11].transAxes)
ax[11].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[11].transAxes)
ax[11].yaxis.set_label_coords(-.15, .5)
#%%
#%%

yr='4x.nc'
yr1='4x'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/nino_temp_control.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/eq_zo_wind_control.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[15].scatter(nino_temp1,nino_Q1)
ax[15].scatter(nino_temp2,nino_Q2)
ax[15].scatter(nino_temp3,nino_Q3)
ax[15].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[15].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[15].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[15].set_xlabel('SST anom. ($^o$C)')
ax[15].set_ylabel('N/m$^2$')
ax[15].set_title('control')
ax[15].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[15].transAxes)
ax[15].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[15].transAxes)
ax[15].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[15].transAxes)

ax[15].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[15].transAxes)
ax[15].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[15].transAxes)
ax[15].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[15].transAxes)
ax[15].yaxis.set_label_coords(-.15, .5)
nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)

from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[16].scatter(nino_temp1,nino_Q1)
ax[16].scatter(nino_temp2,nino_Q2)
ax[16].scatter(nino_temp3,nino_Q3)
ax[16].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[16].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[16].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[16].set_xlabel('SST anom. ($^o$C)')
ax[16].set_ylabel('N/m$^2$')
ax[16].set_title('4xCO$_2$')
ax[16].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[16].transAxes)
ax[16].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[16].transAxes)
ax[16].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[16].transAxes)

ax[16].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[16].transAxes)
ax[16].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[16].transAxes)
ax[16].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[16].transAxes)
ax[16].yaxis.set_label_coords(-.15, .5)

yr='ens1.nc'
yr1='ens1'
yr2='ens1'
yr3='ens1'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens1= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/nino_temp_'+yr2+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens2= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/nino_temp_'+yr3+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens3= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')


nino_temp4=xr.concat([nino_temp_ens1,nino_temp_ens2,nino_temp_ens3],'time')

nino_temp1=nino_temp4.where(nino_temp4>1)
nino_temp2=nino_temp4.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_temp3=nino_temp4.where(nino_temp4<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens1= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/eq_zo_wind_'+yr2+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens2= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC6/eq_zo_wind_'+yr3+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens3= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.concat([nino_Q_ens1,nino_Q_ens2,nino_Q_ens3],'time')


nino_Q1=nino_Q.where(nino_temp4>1)
nino_Q2=nino_Q.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_Q3=nino_Q.where(nino_temp4<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[17].scatter(nino_temp1,nino_Q1)
ax[17].scatter(nino_temp2,nino_Q2)
ax[17].scatter(nino_temp3,nino_Q3)
ax[17].plot(nino_temp4,(nino_temp1*slope1+intercept1),'k')
ax[17].plot(nino_temp4,(nino_temp2*slope2+intercept2),'k')
ax[17].plot(nino_temp4,(nino_temp3*slope3+intercept3),'k')
ax[17].set_xlabel('SST anom. ($^o$C)')
ax[17].set_ylabel('N/m$^2$')
ax[17].set_title('SSP585 (3 ens. mem. )')
ax[17].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[17].transAxes)
ax[17].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[17].transAxes)
ax[17].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[17].transAxes)

ax[17].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[17].transAxes)
ax[17].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[17].transAxes)
ax[17].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[17].transAxes)
ax[17].yaxis.set_label_coords(-.15, .5)
#%%

yr='4x.nc'
yr1='4x'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/IPSL-CM6A-LR/nino_temp_control.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/IPSL-CM6A-LR/eq_zo_wind_control.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[12].scatter(nino_temp1,nino_Q1)
ax[12].scatter(nino_temp2,nino_Q2)
ax[12].scatter(nino_temp3,nino_Q3)
ax[12].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[12].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[12].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[12].set_xlabel('SST anom. ($^o$C)')
ax[12].set_ylabel('N/m$^2$')
ax[12].set_title('control')
ax[12].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[12].transAxes)
ax[12].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[12].transAxes)
ax[12].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[12].transAxes)

ax[12].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[12].transAxes)
ax[12].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[12].transAxes)
ax[12].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[12].transAxes)
ax[12].yaxis.set_label_coords(-.15, .5)


yr='ens1.nc'
yr1='ens1'
yr2='ens1'
yr3='ens1'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/IPSL-CM6A-LR/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens1= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/IPSL-CM6A-LR/nino_temp_'+yr2+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens2= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/IPSL-CM6A-LR/nino_temp_'+yr3+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens3= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')


nino_temp4=xr.concat([nino_temp_ens1,nino_temp_ens2,nino_temp_ens3],'time')

nino_temp1=nino_temp4.where(nino_temp4>1)
nino_temp2=nino_temp4.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_temp3=nino_temp4.where(nino_temp4<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/IPSL-CM6A-LR/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens1= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/IPSL-CM6A-LR/eq_zo_wind_'+yr2+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens2= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/IPSL-CM6A-LR/eq_zo_wind_'+yr3+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens3= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.concat([nino_Q_ens1,nino_Q_ens2,nino_Q_ens3],'time')


nino_Q1=nino_Q.where(nino_temp4>1)
nino_Q2=nino_Q.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_Q3=nino_Q.where(nino_temp4<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[14].scatter(nino_temp1,nino_Q1)
ax[14].scatter(nino_temp2,nino_Q2)
ax[14].scatter(nino_temp3,nino_Q3)
ax[14].plot(nino_temp4,(nino_temp1*slope1+intercept1),'k')
ax[14].plot(nino_temp4,(nino_temp2*slope2+intercept2),'k')
ax[14].plot(nino_temp4,(nino_temp3*slope3+intercept3),'k')
ax[14].set_xlabel('SST anom. ($^o$C)')
ax[14].set_ylabel('N/m$^2$')
ax[14].set_title('SSP585 (3 ens. mem. )')
ax[14].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[14].transAxes)
ax[14].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[14].transAxes)
ax[14].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[14].transAxes)

ax[14].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[14].transAxes)
ax[14].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[14].transAxes)
ax[14].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[14].transAxes)
ax[14].yaxis.set_label_coords(-.15, .5)
#%%

yr='4x.nc'
yr1='4x'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/nino_temp_control.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/eq_zo_wind_control.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[18].scatter(nino_temp1,nino_Q1)
ax[18].scatter(nino_temp2,nino_Q2)
ax[18].scatter(nino_temp3,nino_Q3)
ax[18].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[18].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[18].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[18].set_xlabel('SST anom. ($^o$C)')
ax[18].set_ylabel('N/m$^2$')
ax[18].set_title('control')
ax[18].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[18].transAxes)
ax[18].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[18].transAxes)
ax[18].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[18].transAxes)

ax[18].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[18].transAxes)
ax[18].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[18].transAxes)
ax[18].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[18].transAxes)
ax[18].yaxis.set_label_coords(-.15, .5)
nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')
nino_temp1=nino_temp.where(nino_temp>1)
nino_temp2=nino_temp.where((nino_temp<1) & (nino_temp>-0.75))
nino_temp3=nino_temp.where(nino_temp<-0.75)

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q.to_array()
nino_Q=nino_Q.transpose().squeeze()
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')
nino_Q1=nino_Q.where(nino_temp>1)
nino_Q2=nino_Q.where((nino_temp<1) & (nino_temp>-0.75))
nino_Q3=nino_Q.where(nino_temp<-0.75)

from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[19].scatter(nino_temp1,nino_Q1)
ax[19].scatter(nino_temp2,nino_Q2)
ax[19].scatter(nino_temp3,nino_Q3)
ax[19].plot(nino_temp,(nino_temp1*slope1+intercept1),'k')
ax[19].plot(nino_temp,(nino_temp2*slope2+intercept2),'k')
ax[19].plot(nino_temp,(nino_temp3*slope3+intercept3),'k')
ax[19].set_xlabel('SST anom. ($^o$C)')
ax[19].set_ylabel('N/m$^2$')
ax[19].set_title('4xCO$_2$')
ax[19].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[19].transAxes)
ax[19].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[19].transAxes)
ax[19].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[19].transAxes)

ax[19].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[19].transAxes)
ax[19].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[19].transAxes)
ax[19].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[19].transAxes)
ax[19].yaxis.set_label_coords(-.15, .5)

yr='ens1.nc'
yr1='ens1'
yr2='ens1'
yr3='ens1'

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/nino_temp_'+yr1+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens1= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/nino_temp_'+yr2+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens2= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')

nino_temp=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/nino_temp_'+yr3+'.nc')
nino_temp=nino_temp['thetao']
nino_temp=nino_temp.groupby('time.month')-nino_temp.groupby('time.month').mean('time')
nino_temp=nino_temp.assign_coords(time=range(0,len(nino_temp.time)))
nino_temp_ens3= dsp.bandpass(nino_temp,1/(7*12),1/(1.5*12),dim='time')


nino_temp4=xr.concat([nino_temp_ens1,nino_temp_ens2,nino_temp_ens3],'time')

nino_temp1=nino_temp4.where(nino_temp4>1)
nino_temp2=nino_temp4.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_temp3=nino_temp4.where(nino_temp4<-0.75)


nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/eq_zo_wind_'+yr1+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens1= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/eq_zo_wind_'+yr2+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens2= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.open_dataset('/Users/ullaheede_1/Documents/ENSO_project/MIROC-ES2L/eq_zo_wind_'+yr3+'.nc')
nino_Q=nino_Q['tauu']
nino_Q=nino_Q.groupby('time.month')-nino_Q.groupby('time.month').mean('time')
nino_Q=nino_Q.assign_coords(time=range(0,len(nino_Q.time)))
nino_Q_ens3= dsp.bandpass(nino_Q,1/(7*12),1/(1.5*12),dim='time')

nino_Q=xr.concat([nino_Q_ens1,nino_Q_ens2,nino_Q_ens3],'time')


nino_Q1=nino_Q.where(nino_temp4>1)
nino_Q2=nino_Q.where((nino_temp4<1) & (nino_temp4>-0.75))
nino_Q3=nino_Q.where(nino_temp4<-0.75)


from scipy import stats

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(nino_temp1.dropna("time"),nino_Q1.dropna("time"))
slope2, intercept2, r_value2, p_value, std_err = stats.linregress(nino_temp2.dropna("time"),nino_Q2.dropna("time"))
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(nino_temp3.dropna("time"),nino_Q3.dropna("time"))

ax[20].scatter(nino_temp1,nino_Q1)
ax[20].scatter(nino_temp2,nino_Q2)
ax[20].scatter(nino_temp3,nino_Q3)
ax[20].plot(nino_temp4,(nino_temp1*slope1+intercept1),'k')
ax[20].plot(nino_temp4,(nino_temp2*slope2+intercept2),'k')
ax[20].plot(nino_temp4,(nino_temp3*slope3+intercept3),'k')
ax[20].set_xlabel('SST anom. ($^o$C)')
ax[20].set_ylabel('N/m$^2$')
ax[20].set_title('SSP585 (3 ens. mem. )')
ax[20].text(0.85, 0.95, '$\\mu_a$='+str("%.4f" % slope1), horizontalalignment='center', verticalalignment='center', transform=ax[20].transAxes)
ax[20].text(0.5, 0.95, '$\\mu_a$='+str("%.4f" % slope2), horizontalalignment='center', verticalalignment='center', transform=ax[20].transAxes)
ax[20].text(0.15, 0.95, '$\\mu_a$='+str("%.4f" % slope3), horizontalalignment='center', verticalalignment='center', transform=ax[20].transAxes)

ax[20].text(0.85, 0.1, 'R='+str("%.2f" % r_value1), horizontalalignment='center', verticalalignment='center', transform=ax[20].transAxes)
ax[20].text(0.5, 0.1, 'R='+str("%.2f" % r_value2), horizontalalignment='center', verticalalignment='center', transform=ax[20].transAxes)
ax[20].text(0.15, 0.1, 'R='+str("%.2f" % r_value3), horizontalalignment='center', verticalalignment='center', transform=ax[20].transAxes)
ax[20].yaxis.set_label_coords(-.15, .5)