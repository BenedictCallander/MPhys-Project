import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
import pwlf
from mpl_toolkits.axes_grid1 import make_axes_locatable
#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il

#runtime calculation 
import time


#Computational functions - simultaneous calculations to make use of multi-core CPU
from joblib import Parallel, delayed

#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter

df1 = pd.read_csv("csv/tng33linear.csv")
df2 = pd.read_csv("csv/tng67linear.csv")
df3 = pd.read_csv("csv/tng99linear.csv")

fig,axs = plt.subplots(1,3,figsize=(30,8))

for i in range(3):
    axs[i].set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=15)
    axs[i].set_xlabel("Total Mass [$M_\odot$]",fontsize=15)
    axs[i].grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
    axs[i].tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1,labelsize=15)
    axs[i].set_ylim(-2,2.5)

    
axs[0].set_xlim(9,12)
divider0 = make_axes_locatable(axs[0])  ;cax0 = divider0.append_axes('right', size='5%', pad=0.05)
divider1 = make_axes_locatable(axs[1])  ;cax1 = divider1.append_axes('right', size='5%', pad=0.05)
divider2 = make_axes_locatable(axs[2])  ;cax2 = divider2.append_axes('right', size='5%', pad=0.05)

im0=axs[0].scatter(np.log10(df1['mass']),np.log10(df1['sfr']),c=(df1['slope']),cmap='magma',vmin=-0.02,vmax=0.01)#,vmin=1.5,vmax=8)
im1=axs[1].scatter(np.log10(df2['mass']),np.log10(df2['sfr']),c=(df2['slope']),cmap='magma',vmin=-0.02,vmax=0.01)#,vmin=1.5,vmax=8)
im2=axs[2].scatter(np.log10(df3['mass']),np.log10(df3['sfr']),c=(df3['slope']),cmap='magma',vmin=-0.02,vmax=0.01)#,vmin=1.5,vmax=8)

fig.colorbar(im0, cax=cax0, orientation='vertical')
fig.colorbar(im1, cax=cax1, orientation='vertical')
fig.colorbar(im2, cax=cax2, orientation='vertical')

fig.tight_layout()
fig.subplots_adjust(top=0.89)
fig.suptitle("Redshift progression of galaxy Mass-SFR Relation", fontsize=20)
fig.savefig('slopes_all.png')