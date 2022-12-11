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

def get(path, params = None):
    #utility function for API reading 

    #Make API request - 
    # Path: url to api page 
    #Params - misc ; Headers = api key 
    r = requests.get(path, params=params, headers=headers)

    #HTTP code - raise error if code return is not 200 (success)
    r.raise_for_status()
    
    #detect content type (json or hdf5) - run appropriate download programme
    
    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r


class history:
    def __init__(self,primeID):
        self.primesub = primeID
        
    def gettree(self):
        keepvals = [21,33,50,67,78,91,99]
        fpath = "files/trees/sublink_mpb_{}.hdf5".format(self.primesub)
        with h5py.File(fpath,'r') as f:
            snapshots = list(f['SnapNum'][:])
            subhalos= list (f['SubfindID'][:])
            desclist = [self.primesub] * len(subhalos)
        df = pd.DataFrame({"snapshot":snapshots, "subhalo": subhalos, "primesub": desclist})
        df = df[df['snapshot'].isin(keepvals)]
        self.df = df
        return df
    def get_globals(self):
        df = self.df
        snaps = list(df['snapshot'])
        subs = list(df['subhalo'])
        mass = [];sfr=[];met = []
        for i in range(len(snaps)):
            snap = snaps[i]
            sub = subs[i]
            url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/".format(snap,sub)
            subhalo = get(url)
            mass.append(subhalo['mass_log_msun'])
            met.append(subhalo['gasmetallicity'])
            sfr.append(subhalo['sfr'])
        df = pd.DataFrame({"snapshot":snaps, "subhalo":subs,"mass":mass,"sfr":sfr,"met":met})
        fpath = "files/historycutouts/evdir_{}/historydata_{}.csv".format(self.primesub,self.primesub)
        df.to_csv()
            
        
    
df = pd.read_csv("traceids.csv")
ids = list(df['id'])

def dofunc(i):
    obj = history(i)
    df = obj.gettree()
    obj.get_globals()
    print("done for {}".format(i))

returns = Parallel(n_jobs=20)(delayed(dofunc)(i) for i in ids)
# point of no return was 572840
#
