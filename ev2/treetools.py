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
        df.to_csv(fpath)
        self.maindf = df
        return df
    def MSPLOT(self):
        #snapshots = [1,5,6,10,21,33,40,50,52,55,67,78,91,99]
        #linesnaps = [1,5,6,10,21,33,40,50,52,55,67,78,91]
        snapshots = [21,33,50,67,78,91,99]
        linesnaps = [21,33,50,67,78,91,99]
        linedf = self.maindf[self.maindf['snapshot'].isin(linesnaps)]
        self.maindf = self.maindf[self.maindf['snapshot'].isin(snapshots)]
        #print(self.maindf)
        basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'
        fields = ['SubhaloMass','SubhaloSFRinRad','SubhaloSFR','SubhaloGasMetallicitySfrWeighted']
        subhalos99 = il.groupcat.loadSubhalos(basePath,33,fields=fields)
        mass =subhalos99['SubhaloMass'] * 1e10 / 0.704
        sfr = list((subhalos99['SubhaloSFR']))
        plt.figure(figsize=(20,12))
        plt.yscale('log')
        plt.title("evolution of largest progenitor to subhalo {} at z=0: ".format(self.primesub))
        plt.plot(np.log10(mass),sfr, 'g+',label = 'Snapshot 99 subhalos', zorder=1)
        plt.plot(linedf['mass'],linedf['sfr'], 'r--',label = 'path of progenitors to subhalo {}'.format(self.primesub),zorder=2)
        plt.scatter(self.maindf['mass'],self.maindf['sfr'],c = self.maindf['snapshot'], cmap = 'magma',vmin=1,vmax=99,zorder=3)
        plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
        plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.colorbar(label = 'Snapshot')
        plt.ylabel('log(SFR)')
        plt.xlabel('Mass (log10 Msun)')
        plt.legend(loc='upper right')
        #'MSPLOT/MS_evolution_{}.png'.format(self.subID)
        plt.savefig('msev_{}.png'.format(self.primesub))
        plt.close()
        
    
df = pd.read_csv("traceids.csv")
ids = list(df['id'])

def dofunc(i):
    try:
        obj = history(i)
        df = obj.gettree()
        obj.get_globals()
        print("done for {}".format(i))
    except FileNotFoundError as e:
        print(e)
    except IndexError as e:
        print(e)

'''
plt.figure(figsize=(20,12))
plt.title("snapshot evolution of properties")
plt.plot(df['snapshot'],df['sfr'],'r-')
plt.xlabel("Snapshot Number")
plt.ylabel("sfr")
plt.savefig("sfrev.png")
plt.close()
'''
returns = Parallel(n_jobs=20)(delayed(dofunc)(i) for i in ids)
# point of no return was 572840
#
