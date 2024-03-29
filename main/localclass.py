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
met33 = list(subhalos33['SubhaloGasMetallicitySfrWeighted']); met67 = list(subhalos67['SubhaloGasMetallicitySfrWeighted'])
met99 = list(subhalos99['SubhaloGasMetallicitySfrWeighted'])

print(np.mean(masses33))
print(np.mean(masses67))
print(np.mean(masses99))

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
divider0 = make_axes_locatable(ax0); divider1 = make_axes_locatable(ax1);divider2 = make_axes_locatable(ax2)
cax0 = divider0.append_axes('right', size='5%', pad=0.05)
cax1 = divider1.append_axes('right', size='5%', pad=0.05)
cax2 = divider2.append_axes('right', size='5%', pad=0.05)


plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

fig,axs = plt.subplots(nrows = 1, ncols = 3,figsize=(30,8))

#set scales to log
#y
axs[0].set_yscale('log');axs[1].set_yscale('log');axs[2].set_yscale('log')
#x

axs[0].set_title("Snapshot 33: z=2",fontsize=20)
axs[0].plot(np.log10(masses33),sfr33, 'g+')
axs[0].set_xlabel("Total Mass [$M_\odot$]",fontsize=20);axs[0].set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
for i in range(3):
    axs[i].grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
    axs[i].tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)




axs[1].set_title("Snapshot 67: z=0.5",fontsize=20)
axs[1].plot(np.log10(masses67),sfr67, 'g+')
axs[1].set_xlabel("Total Mass [$M_\odot$]",fontsize=20);axs[1].set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
#axs[1].grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
#axs[1].tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)

axs[2].set_title("Snapshot 99: z=0",fontsize=20)
axs[2].plot(np.log10(masses99),sfr99, 'g+')
axs[2].set_xlabel("Total Mass [$M_\odot$]",fontsize=20);axs[2].set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
#axs[2].grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
#axs[2].tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)

fig.tight_layout()
fig.subplots_adjust(top=0.89)
fig.suptitle("Redshift progression of galaxy Mass-SFR Relation", fontsize=20)
fig.savefig('localclass_all.png')

'''
