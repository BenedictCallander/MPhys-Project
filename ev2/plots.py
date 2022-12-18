import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
import pwlf

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



df = pd.read_csv("tng67KPCslopeslin.csv")

    
    
print(min(df['slope']))
print(max(df['slope']))
print(0.5*np.mean(df['slope']))
print(1.5*np.mean(df['slope']))

plt.figure(figsize=(20,12))

plt.yscale('log')

plt.title("Metallicity Gradients of TNG50-1 Main Sequence Subhalos at z=0.5")
plt.xlabel("Total Mass [$M_\odot$]",fontsize=20)
plt.ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)

plt.scatter(df['mass'],(df['sfr']),c=(df['slope']),cmap='magma',vmin=-0.5,vmax=0.3)
plt.colorbar(label="Metallicity Gradient [dex/Kpc]")

plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.savefig("slopechar67.png")
plt.close()