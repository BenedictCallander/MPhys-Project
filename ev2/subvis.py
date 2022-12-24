import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
import pwlf
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as pat


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


plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

dens = pd.read_csv("dens.csv")
met = pd.read_csv("met.csv")
sfr = pd.read_csv("sfr.csv")
hist = pd.read_csv("historgram.csv")
plt.style.use('dark_background')
fig,((ax0,ax1),(ax2,ax3)) = plt.subplots(nrows = 2, ncols = 2,figsize=(20,16))
divider0 = make_axes_locatable(ax0)  ;cax0 = divider0.append_axes('right', size='5%', pad=0.05)
divider1 = make_axes_locatable(ax1)  ;cax1 = divider1.append_axes('right', size='5%', pad=0.05)
divider2 = make_axes_locatable(ax2)  ;cax2 = divider2.append_axes('right', size='5%', pad=0.05)
#divider3 = make_axes_locatable(ax3)  ;cax3 = divider2.append_axes('right', size='5%', pad=0.05)

ax0.set_xlabel('$\Delta x$ [kpc/h]')
ax0.set_ylabel('$\Delta y$ [kpc/h]')


annulus0 =pat.Circle((0, 0), 30, color='g', fill=False,lw=5)
ax0.set_title("Subhalo Gas Density")
im0=ax0.scatter(dens['x'],dens['y'],c=(np.log10(dens['m'])),cmap='inferno', vmin=(min(np.log10(dens['m']))),vmax =(0.7*max(np.log10(dens['m']))))
ax0.add_patch(annulus0)

annulus1 =pat.Circle((0, 0), 30, color='g', fill=False,lw=5)
ax1.set_title("Subhalo Surface Metallicity")
im1=ax1.scatter(met['x'],met['y'],c=(np.log10(met['met'])),cmap='inferno', vmin=(min(np.log10(met['met']))),vmax =(0.7*max(np.log10(met['met']))))
ax1.add_patch(annulus1)

annulus2 =pat.Circle((0, 0), 30, color='g', fill=False,lw=5)
ax2.set_title("Subhalo Surface SFR")
im2=ax2.scatter(sfr['x'],sfr['y'],c=(np.log10(sfr['sfr'])),cmap='inferno')# vmin=(min(np.log10(df_valid['m']))),vmax =(0.7*max(np.log10(df_valid['m']))))
ax2.add_patch(annulus2)

ax3.set_title("Subhalo Metallicity Profile")
im3=ax3.hist2d(hist['rad'],12+np.log10(hist['met']),bins=[200,200], weights=1/hist['sfr'],cmap='PuOr')
ax3.set_xlabel("Radius (Normalised Code Units)")
ax3.set_ylabel("12+$log_{10}$ $(O/H)$")





fig.colorbar(im0, cax=cax0, orientation='vertical')
fig.colorbar(im1, cax=cax1, orientation='vertical')
fig.colorbar(im2, cax=cax2, orientation='vertical')

fig.tight_layout()

fig.savefig("subvis.png")   
