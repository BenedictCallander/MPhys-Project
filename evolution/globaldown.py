#
# Load progenitor subhalos to subset of TNG99 subhalos and compute their metallicity parameters
#
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
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter
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
    def __init__(self, descendant):
        self.startsub = descendant
        
        self.mpb = "files/binary/trees/99trees/sublink_mpb_{}.hdf5".format(self.startsub)
        keepvals = [21,33,50,67,78,91,99]
        self.length = keepvals
        with h5py.File(self.mpb,'r') as f:
            snapshots = list(f['SnapNum'][:])
            subhalos= list (f['SubfindID'][:])
        df = pd.DataFrame({
            "snapshots": snapshots,
            "subhalos":subhalos
        })
        
        df = df[df['snapshots'].isin(keepvals)]
        df.to_csv("files/historycutouts/evdir_{}/treedata_{}.csv".format(self.startsub,self.startsub))
        self.target_snaps= list(df['snapshots'])
        self.target_subhalos = list(df['subhalos'])
        #print(df)
    
    def cutoutdownload(self):
        mass = []; sfr = []; met = []; Rhalf = []; stellarphotometrics = []
        snaps = []; subs = []
        for i in range(6):
            snap = self.target_snaps[i]
            snaps.append(snap)
            sub = self.target_subhalos[i]
            subs.append(sub)
            url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/".format(snap,sub)
            subhalo = get(url)
            mass.append(subhalo['mass_log_msun'])
            sfr.append(subhalo['sfr'])
            met.append(subhalo['gasmetallicity'])
            Rhalf.append(subhalo['halfmassrad'])
            stellarphotometrics.append(subhalo['stellarphotometricsrad'])
        fpath = "files/historycutouts/evdir_{}/subdata_{}.hdf5".format(self.startsub,self.startsub)
        with h5py.File(fpath, 'w') as f:
            # Write the lists to the HDF5 file
            f['list1'] = mass; f['list2'] = sfr ; f['list3'] = met ; f['list4'] = Rhalf; f['list5'] = stellarphotometrics 
            f['list6'] = snaps ; f['list7'] = subs
            # Set the keywords for the lists
            f['list1'].attrs['keyword'] = 'mass'
            f['list2'].attrs['keyword'] = 'sfr'
            f['list3'].attrs['keyword'] = 'met'
            f['list4'].attrs['keyword'] = 'Rhalf'
            f['list5'].attrs['keyword'] = 'stellarph'
            f['list6'].attrs['keyword'] = 'snaps'
            f['list7'].attrs['keyword'] = 'subs'
            
df = pd.read_csv("traceids.csv")
ids = list(df['id'])               
            
def down(i):
    sub = history(i)
    return print("subahlo {} done".format(i))

returns = Parallel(n_jobs= 20)(delayed(down)(j) for j in ids)