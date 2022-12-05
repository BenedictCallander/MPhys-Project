# Subhalo.py 
# \-> contains subhalo class and subsequent analysis function (which runs all desired operations on each subhalo object)
# Created:17/11/2022 
# Author: Benedict Callander 
# Respository https://github.com/btcallander/MPhys-Project (private)
#

#Plotting, numerical functions and dataset manipulation
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import pwlf

#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il

#Own module containing utility functions 
import BCUTILS

#runtime calculation 
import time

#Computational functions - simultaneous calculations to make use of multi-core CPU
from joblib import Parallel, delayed

#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter

Msun = 1.98847e30
#set basic constants during initialisation for easy 
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
pd.options.mode.chained_assignment = None  # default='warn'

basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'
fields = ['SubhaloMass','SubhaloSFRinRad','SubhaloSFR','SubhaloGasMetallicitySfrWeighted']
subhalos33 = il.groupcat.loadSubhalos(basePath,33,fields=fields)
subhalos67 = il.groupcat.loadSubhalos(basePath,67,fields=fields)
subhalos99 = il.groupcat.loadSubhalos(basePath,99,fields=fields)


masses33 =subhalos33['SubhaloMass'] * 1e10 / 0.704 ;masses67 =subhalos67['SubhaloMass'] * 1e10 / 0.704 
masses99 =subhalos99['SubhaloMass'] * 1e10 / 0.704

sfr33 = list((subhalos33['SubhaloSFR'])); sfr67 = list((subhalos67['SubhaloSFR'])) ;sfr99 = list((subhalos99['SubhaloSFR']))

'''
plt.figure(figsize=(20,12))
plt.plot((masses99), sfr99, 'g+')
plt.yscale('log')
plt.xscale('log') 
plt.xlabel("Total Mass [$M_\odot$]",fontsize=20);plt.ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.savefig('m_sfr_99_local.png')
plt.close()
'''

plt.figure(figsize=(20,12))
plt.plot( sfr99,(masses99), 'g+')
plt.yscale('log')
plt.xscale('log') 
plt.ylabel("Total Mass [$M_\odot$]",fontsize=20);plt.xlabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.savefig('sfr_m_99_local.png')
plt.close()



'''
divider0 = make_axes_locatable(ax0); divider1 = make_axes_locatable(ax1);divider2 = make_axes_locatable(ax2)
cax0 = divider0.append_axes('right', size='5%', pad=0.05)
cax1 = divider1.append_axes('right', size='5%', pad=0.05)
cax2 = divider2.append_axes('right', size='5%', pad=0.05)


plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

fig,(ax0,ax1,ax2) = plt.subplots(nrows = 1, ncols = 3,figsize=(30,8))

#set scales to log
#y
ax0.set_yscale('log');ax1.set_yscale('log');ax2.set_yscale('log')
#x
ax0.set_xscale('log');ax1.set_xscale('log');ax2.set_xscale('log')

ax0.set_title("Snapshot 33: z=2",fontsize=20)
ax0.plot(masses33,sfr33, 'g+')
ax0.set_xlabel("Total Mass [$M_\odot$]",fontsize=12);ax0.set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=12)
ax0.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
ax1.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
ax2.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)

ax0.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
ax1.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
ax2.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)



ax1.set_title("Snapshot 67: z=0.5",fontsize=20)
ax1.plot(masses67,sfr67, 'g+')
ax1.set_xlabel("Total Mass [$M_\odot$]",fontsize=12);ax0.set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=12)
ax1.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
ax1.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)

ax2.set_title("Snapshot 99: z=0",fontsize=20)
ax2.plot(masses99,sfr99, 'g+')
ax2.set_xlabel("Total Mass [$M_\odot$]",fontsize=12);ax0.set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=12)
ax2.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
ax2.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)

fig.tight_layout()
fig.subplots_adjust(top=0.89)
fig.suptitle("Redshift progression of galaxy Mass-SFR Relation", fontsize=20)
fig.savefig('localclass_all.png')

'''
