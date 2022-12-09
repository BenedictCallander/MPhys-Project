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
    def __init__(self,final):
        self.primesub = final
        
    def gettree(self):
        locals = [33,67,99]
        fpath = "files/binary/trees/99trees/sublink_mpb_{}.hdf5".format(self.primesub)
        with h5py.File (fpath,'r') as f:
            snapshots = list(f['SnapNum'][:])
            subhalos= list (f['SubfindID'][:])
            
        df = pd.DataFrame({"snapshots": snapshots, "subhalos":subhalos })
        df = df[df['snapshots'].isin(locals)]
        
        writepath = "files/historycutouts/evdir_{}/locals_{}.csv".format(self.primesub,self.primesub)
        df.to_csv(writepath)
        return print("done for {}".format(self.primesub))


df = pd.read_csv("traceids.csv")
ids = list(df['id'])
def writesubs(i):
    filedo = history(i)
    filedo.gettree()


    

returns = Parallel(n_jobs=20)(delayed(writesubs)(i)for i in ids)

    