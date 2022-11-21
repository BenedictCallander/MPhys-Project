import os 
import time  # runtime calculation import numpy as np #data handling
from random import random
from re import sub  # http logging for debugging purpouses
import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests  # obtain data from API server
import h5py
from joblib import Parallel,delayed

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

url96762 = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/96762/"
sub = get(url96762)

#mpb1 - main progenitor branch 1 
mpb1 ="hdf5/sublink_mpb_0.hdf5"
#get(sub['trees']['sublink_mpb'] )
 # file saved, mpb1 contains the filename

f = h5py.File(mpb1,'r')
#print(f.keys())
#print(f['SnapNum'][:])
mpb2 ="hdf5/lhalotree_mpb_0.hdf5"
#get( sub['trees']['lhalotree_mpb'] )
# # file saved, mpb2 contains the filename

f = h5py.File(mpb2, 'r')
#print(len(f['SnapNum'][:]))

with h5py.File(mpb2,'r') as f:
    pos = f['SubhaloPos'][:]
    snapnum = f['SnapNum'][:]
    subid = f['SubhaloNumber'][:]
    GM = f['SubhaloGasMetallicity'][:]
    SFR = f['SubhaloSFR'][:]
    SM = f['SubhaloStarMetallicity'][:]

snapnum = list(snapnum); subid = list(subid)
snapnum.reverse(); subid.reverse()
print(subid); print(snapnum)

sub96762_ids = pd.DataFrame({
    "snapnum": snapnum,
    "subID": subid
})



def url_constructor(snap,subid):
    baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/".format(snap,subid)
    return baseurl

def sub_trace(i):
    url = url_constructor(snapnum[i], subid[i])
    sub = get(url)
    met = sub['gasmetallicity']+sub['starmetallicity']
    mass = sub['mass_log_msun']
    sfr = sub['sfr']
    print("snap{} done!".format(snapnum[i]))
    return (met,mass,sfr,subid[i],snapnum[i])
returns = Parallel(n_jobs= 15)(delayed(sub_trace)(i) for i in range(1,99))
df=pd.DataFrame(returns,columns=['met','mass','sfr','id','snapshot'])
df.to_csv("sublinks_0.csv")

plt.figure(figsize = (20,12))
plt.plot(df['snapshot'], df['mass'], 'g-')
plt.plot(df['snapshot'], df['sfr'],'r-')
plt.xlabel("snapshot", fontsize = 20)
plt.ylabel("Mass (Log_MSUN)", fontsize = 20)
plt.savefig("sub2.png")
plt.close()

df2 = pd.read_csv("csv/tng99subhalos.csv")
sample = df2.sample(frac=0.1, replace=False)


plt.figure(figsize = (20,12))
plt.plot(df2['mass'], df2['sfr'], 'g.', ms = 3)
plt.scatter(df['mass'], df['sfr'],c = df['snapshot'], cmap = 'viridis', s=20)
plt.yscale('log')
plt.xlabel("snapshot", fontsize = 20)
plt.ylabel("Mass (Log_MSUN)", fontsize = 20)
plt.savefig("sub4.png")
plt.close()