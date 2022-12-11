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


class history:
    def __init__(self,primeID):
        self.primesub = primeID
        
    def gettree(self):
        keepvals = [21,33,50,67,78,91,99]
        fpath = "/home/AstroPhysics-Shared/PROJECTS/MPhys_Schady/Projects/TNGmetgrads/BC763/MPhys-Project/evolution/files/binary/trees/99trees/sublink_mpb_{}.hdf5".format(self.primesub)
        with h5py.File(fpath,'r') as f:
            snapshots = list(f['SnapNum'][:])
            subhalos= list (f['SubfindID'][:])
            desclist = [self.primesub] * len(subhalos)
        df = pd.DataFrame({"snapshot":snapshots, "subhalo": subhalos, "primesub": desclist})
        df = df[df['snapshot'].isin(keepvals)]
        return df
    
df = pd.read_csv("traceids.csv")
ids = list(df['id'])

dataframes = []
for i in ids:
    obj = history(i)
    df = obj.gettree()
    dataframes.append(df)
    print("done for {}".format(i))
combined = pd.concat(dataframes)
combined.to_csv("all.csv")


#
# point of no return was 572840
#
