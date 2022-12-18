import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
import pwlf
from mpl_toolkits.axes_grid1 import make_axes_locatable


#hdf5 binary file manipulation
import h5py

#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il

#Own module containing utility functions 

#runtime calculation 
import time

#Computational functions - simultaneous calculations to make use of multi-core CPU
from joblib import Parallel, delayed

#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()


df1 = pd.read_csv("csv/tng33KPCslopesboth.csv")
df2= pd.read_csv("csv/tng67KPCslopesboth.csv")
df3 = pd.read_csv("csv/tng99KPCslopesboth.csv")
'''

df1= pd.read_csv("tng33MS2slopes.csv")
df2= pd.read_csv("tng67MS2slopes.csv")
df3= pd.read_csv("tng99MS2slopes.csv")
'''

print("\n33\n")
print(min(df1['slope1']))
print(max(df1['slope1']))
print(np.mean(df1['slope1']))
print("\n67\n")
print(min(df2['slope1']))
print(max(df2['slope1']))
print(np.mean(df2['slope1']))
print("\n99\n")
print(min(df3['slope1']))
print(max(df3['slope1']))
print(np.mean(df3['slope1']))



plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

fig,axs = plt.subplots(nrows = 1, ncols = 3,figsize=(30,8))

divider0 = make_axes_locatable(axs[0])  ;cax0 = divider0.append_axes('right', size='5%', pad=0.05)
divider1 = make_axes_locatable(axs[1])  ;cax1 = divider1.append_axes('right', size='5%', pad=0.05)
divider2 = make_axes_locatable(axs[2])  ;cax2 = divider2.append_axes('right', size='5%', pad=0.05)

for i in range(3):
    axs[i].grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
    axs[i].tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
    axs[i].set_yscale('log')
    axs[i].set_xlabel("Total Mass [$M_\odot$]",fontsize=20)
    axs[i].set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)

axs[0].set_title("Snapshot 33: z=2",fontsize=20)
axs[1].set_title("Snapshot 67: z=0.5",fontsize=20)
axs[2].set_title("Snapshot 99: z=0",fontsize=20)

im0=axs[0].scatter(df1['mass'],(df1['sfr']),c=(df1['slope1']),cmap='magma',vmin=-0.6,vmax=0.2)
im1=axs[1].scatter(df2['mass'],(df2['sfr']),c=(df2['slope1']),cmap='magma',vmin=-0.6,vmax=0.2)
im2=axs[2].scatter(df3['mass'],(df3['sfr']),c=(df3['slope1']),cmap='magma',vmin=-0.6,vmax=0.2)

fig.colorbar(im0, cax=cax0, orientation='vertical')
fig.colorbar(im1, cax=cax1, orientation='vertical')
fig.colorbar(im2, cax=cax2, orientation='vertical')

fig.tight_layout()
fig.subplots_adjust(top=0.89)
fig.suptitle("Redshift progression of galaxy Mass-SFR Relation", fontsize=20)
fig.savefig('inner_slopes_all.png')