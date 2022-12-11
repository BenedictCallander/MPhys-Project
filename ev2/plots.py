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



df = pd.read_csv("alldatadone.csv")
df.dropna()
slopes1 =list(df['slope1'])
print(np.mean(slopes1))
slopes2 =list(df['slope2'])
print(np.mean(slopes2))
slopes3 = list(df['slope3'])
slopechar1=[]
slopechar2=[]
for i in range(len(slopes1)):
    slopechar1.append(slopes2[i]-slopes1[i])
    slopechar2.append(slopes3[i]-slopes2[i])
    
    
    

df = df[df['snapshot'].isin([78])]
plt.figure(figsize=(20,12))
plt.scatter(df['mass'],(df['sfr']),c=(df['slope2']-df['slope1']) ,cmap = 'magma')
plt.yscale('log')
plt.xlabel("Log10M")
plt.colorbar(label="Outerslope")
plt.ylabel("12+log10(O/H)")
plt.savefig("slopechar1.png")
plt.close()