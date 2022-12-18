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
#subhalos33 = il.groupcat.loadSubhalos(basePath,33,fields=fields)
#subhalos67 = il.groupcat.loadSubhalos(basePath,67,fields=fields)
subhalos99 = il.groupcat.loadSubhalos(basePath,99,fields=fields)

df = pd.read_csv('csv/tng99slopes.csv')
slopes= list(df['slope'])

#masses33 =subhalos33['SubhaloMass'] * 1e10 / 0.704 ;masses67 =subhalos67['SubhaloMass'] * 1e10 / 0.704 
masses99 =subhalos99['SubhaloMass'] * 1e10 / 0.704

#sfr33 = list((subhalos33['SubhaloSFR'])); sfr67 = list((subhalos67['SubhaloSFR'])) 
sfr99 = list((subhalos99['SubhaloSFR']))
#met33 = list(subhalos33['SubhaloGasMetallicitySfrWeighted']); met67 = list(subhalos67['SubhaloGasMetallicitySfrWeighted'])
met99 = list(subhalos99['SubhaloGasMetallicitySfrWeighted'])
df1 = pd.DataFrame({'mass':masses99, 'sfr':sfr99,'met':met99})
df1 = df[df['sfr']>0]
print(len(df1['sfr']))

plt.figure(figsize=(20,12))
plt.scatter(np.log10(df1['mass']),12+np.log10(df1['met']), c = slopes, cmap = 'magma',vmin = -0.5,vmax=0.5)
#plt.xscale('log') 
plt.xlabel("Total Mass [$M_\odot$]",fontsize=20);plt.ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.savefig('MZR99slope.png')
plt.close()