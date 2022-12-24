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



subhalos99 = il.groupcat.loadSubhalos(basePath,99,fields=fields)
masses99 =subhalos99['SubhaloMass'] * 1e10 / 0.704
sfr99 = list((subhalos99['SubhaloSFR']))
met99 = list(subhalos99['SubhaloGasMetallicitySfrWeighted'])

subhalos33 = il.groupcat.loadSubhalos(basePath,33,fields=fields)
masses33 =subhalos33['SubhaloMass'] * 1e10 / 0.704 
sfr33 = list((subhalos33['SubhaloSFR']))
met33 = list(subhalos33['SubhaloGasMetallicitySfrWeighted'])

subhalos67 = il.groupcat.loadSubhalos(basePath,67,fields=fields)
masses67 =subhalos67['SubhaloMass'] * 1e10 / 0.704 
sfr67 = list((subhalos67['SubhaloSFR'])) 
met67 = list(subhalos67['SubhaloGasMetallicitySfrWeighted'])

df33 = pd.DataFrame({'mass':masses33,'sfr':sfr33,'met':met33})
df33 = df33[df33['sfr']>0]

df67 = pd.DataFrame({'mass':masses67,'sfr':sfr67,'met':met67})
df67 = df67[df67['sfr']>0]

df99 = pd.DataFrame({'mass':masses99,'sfr':sfr99,'met':met99})
df99 = df99[df99['sfr']>0]


'''

plt.figure(figsize=(20,12))
plt.title("Fundamental Metallicity Relation for TNG50-1 Snapshot 67",fontsize=25)
plt.plot(np.log10(df['mass']),12+np.log10(df['met']), 'g+')
#plt.scatter(df['mass'], 12+np.log10(df['met']),c=np.log10(df['sfr']),cmap='magma')
#plt.xscale('log')
plt.xlabel("Subhalo Mass [$M_\odot$]",fontsize=20);plt.ylabel(r"12+$log_{10}$ $(\frac{O}{H})$ : [SFR weighted]",fontsize=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#plt.colorbar(label = "$Log_{10}(SFR)$")
plt.tight_layout()
plt.savefig('MZR67.png')
plt.close()
'''


fig,axs = plt.subplots(nrows = 1, ncols = 3,figsize=(30,8),dpi=500)



for i in range(3):
    axs[i].grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
    axs[i].tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
    axs[i].set_xscale('log')
    axs[i].set_yscale('log')
    axs[i].set_xlabel("Total Mass [$M_\odot$]",fontsize=20)
    axs[i].set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)

axs[0].set_title("Snapshot 33: z~2",fontsize=20)
axs[1].set_title("Snapshot 67: z~0.5",fontsize=20)
axs[2].set_title("Snapshot 99: z=0",fontsize=20)

im0=axs[0].plot(df33['mass'],df33['sfr'],'g+')
im1=axs[1].plot(df67['mass'],df67['sfr'],'g+')
im2=axs[2].plot(df99['mass'],df99['sfr'],'g+')



fig.tight_layout()
fig.subplots_adjust(top=0.89)
fig.suptitle("TNG50-1 Subhalo Population ", fontsize=20)
fig.savefig('aa.png')


end = time.time()
print("complete time {}".format(end-start))
'''
divider0 = make_axes_locatable(axs[0])  ;cax0 = divider0.append_axes('right', size='5%', pad=0.05)
divider1 = make_axes_locatable(axs[1])  ;cax1 = divider1.append_axes('right', size='5%', pad=0.05)
divider2 = make_axes_locatable(axs[2])  ;cax2 = divider2.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im0, cax=cax0, orientation='vertical')
fig.colorbar(im1, cax=cax1, orientation='vertical')
fig.colorbar(im2, cax=cax2, orientation='vertical')
'''