import os 
import time  # runtime calculation import numpy as np #data handling
from random import random
from re import sub  # http logging for debugging purpouses
import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests  # obtain data from API server
import pwlf
import h5py
from joblib import Parallel, delayed
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter


headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
pd.options.mode.chained_assignment = None  # default='warn'

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


class subtree:
    def __init__(self, knownID, knownSnap):
        self.knownID = knownID
        self.knownSnap = knownSnap
        
        self.baseurl = "https://www.tng-project.org/api"
        self.known_url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}".format(str(self.knownSnap), str(self.knownID))
    
    def get_trees(self):
        sub = get(self.known_url) 
        #mpb1 = get(sub['trees']['sublink_mpb'])
        mpb2 = "lhalotree_mpb_11.hdf5"
        #get( sub['trees']['lhalotree_mpb'])
        #self.mpb1 = mpb1
        self.mpb2 = mpb2
        
    def get_IDS(self, filename):
        
        f = h5py.File(self.mpb2, 'r')
        with h5py.File(self.mpb2,'r') as f:
            snapnum = f['SnapNum'][:]
            subid = f['SubhaloNumber'][:]
        snapnum = list(snapnum); subid = list(subid)
        snapnum.reverse(); subid.reverse()
        df = pd.DataFrame({
            'snapshot': snapnum,
            'id': subid
        })
        df.to_csv(filename)
        self.histdf = df
        self.snapnum=snapnum
        self.subid= subid
    def get_history(self,i):
        url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/".format(self.snapnum[i], self.subid[i])
        sub = get(url)
        met = sub['gasmetallicity']+sub['starmetallicity']
        mass = sub['mass_log_msun']
        sfr = sub['sfr']
        print("snap{} done!".format(self.snapnum[i]))
        return (met,mass,sfr,self.subid[i],self.snapnum[i])

    def history_plot(self, file, property):
        df = pd.read_csv(file)
        plt.figure(figsize=(20,12))
        
        plt.title("Redshift evolution for {}: subhalo {}".format(str(property),str(self.knownID)))
        if property is "mass":
            plt.plot(df['snapshot'], df['mass'],'r-')
            plt.ylabel("Mass ($log_{10} M_{sun}$)")
        elif property is "sfr":
            plt.yscale('log')
            plt.plot(df['snapshot'], df['sfr'],'b-')
            plt.ylabel("$log_10$(SFR)")
        elif property is "metallicity":
            plt.plot(df['snapshot'], 12+np.log10(df['met']),'r-')
            plt.ylabel("$12 + log_{10}(O/H)$")
        plt.xlabel("Snapshot Number")
        pngname= "{}_ev".format(property)
        plt.savefig(pngname)
        plt.close()
            
        
        
        
        
        
        
subtr = subtree(11,99)
subtr.get_trees()
subtr.get_IDS("sub11trees.csv")

returns = Parallel(n_jobs= 15)(delayed(subtr.get_history)(i) for i in range(1,99))
df=pd.DataFrame(returns,columns=['met','mass','sfr','id','snapshot'])
file = "sublinks_11.csv"
df.to_csv(file)

subtr.history_plot(file,"mass")
subtr.history_plot(file,"metallicity")
subtr.history_plot(file,"sfr")
