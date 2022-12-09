# Metallicityev.py 
# \-> script containing Classes subsequent functions to study the Metallicity evolution of a subhalo's metallicity gradient through the IllustrisTNG snapshots 
# Created:17/11/2022 
# Author: Benedict Callander 
# Respository https://github.com/btcallander/MPhys-Project (private)
#


#Plotting, numerical functions and dataset manipulation
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
import pwlf
import os 
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

def varget(path,dir, params = None):
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
        filename = dir + r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

'''
Input parameters

File containing IDS of all TNG99 subhalos to trace metallicity evolution

'''
df = pd.read_csv("traceids2.csv")
ids = list(df['id'])
'''
for i in ids:
    dir_name = "files/binary/historycutouts/evdir_{}".format(i)
    os.makedirs(dir_name)
'''

class history:
    def __init__(self, descendant):
        self.startsub = descendant
        
        self.mpb = "files/binary/trees/99trees/sublink_mpb_{}.hdf5".format(self.startsub)
        keepvals = [21,50,78,91]
        self.length = keepvals
        with h5py.File(self.mpb,'r') as f:
            snapshots = list(f['SnapNum'][:])
            subhalos= list (f['SubfindID'][:])
        df = pd.DataFrame({
            "snapshots": snapshots,
            "subhalos":subhalos
        })
        
        df = df[df['snapshots'].isin(keepvals)]
        self.target_snaps= list(df['snapshots'])
        self.target_subhalos = list(df['subhalos'])
        #print(df)
    
    def cutoutdownload(self):
        for i in range(4):
            snap = self.target_snaps[i]
            sub = self.target_subhalos[i]
            url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/".format(snap,sub)
            temp = get(url)
            cutout_request = {'gas':'Coordinates,Masses,GFM_Metallicity,StarFormationRate,Velocities'}
            cutout = varget(temp['cutouts']['subhalo'],"files/binary/cutouts/".format(self.startsub),cutout_request)
            
def down(i):
    sub = history(i)
    sub.cutoutdownload()
    return print("subahlo {}".format(i))

returns = Parallel(n_jobs= 20)(delayed(down)(j) for j in ids)
