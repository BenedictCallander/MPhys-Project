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
#import pwlf

#hdf5 binary file manipulation
import h5py

#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il

#Own module containing utility functions 
import BCUTILS

#runtime calculation 
import time

#Computational functions - simultaneous calculations to make use of multi-core CPU
from joblib import Parallel, delayed

#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
pd.options.mode.chained_assignment = None  # default='warn'

class UTILITY:
    def get(path, params = None):
        r'''
        func get(path, params = None)
        
        Utility function to read data from API using http requests: custom modification expanding upon requests.get():
        
        
        Expansions:
        
        Content Filtering: built to read information from illustris TNG API:
        Error Raising: include raise_for_status function to display all error codes (expect HTTP return 200 to indicate all OK)
        
        Valid data types =
        
        ->application/json
        
        -> .hdf5 
        
        
        '''
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
            filename = "hdf5/"+r.headers['content-disposition'].split("filename=")[1]
            with open(filename, 'wb') as f:
                f.write(r.content)
            return filename # return the filename string

        return r
    
    def treeget(path, params = None):
        r'''
        func get(path, params = None)
        
        Utility function to read data from API using http requests: custom modification expanding upon requests.get():
        
        
        Expansions:
        
        Content Filtering: built to read information from illustris TNG API:
        Error Raising: include raise_for_status function to display all error codes (expect HTTP return 200 to indicate all OK)
        
        Valid data types =
        
        ->application/json
        
        -> .hdf5 
        
        
        '''
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
            filename = "trees/"+r.headers['content-disposition'].split("filename=")[1]
            with open(filename, 'wb') as f:
                f.write(r.content)
            return filename # return the filename string

        return r

    def line(m,x,b):
        '''
        straight line function y=m*x+b
        '''
        y = 10**((m*x)+b)
        return y 
    
    def linear_fit(a,x,b):
        f = (a*x)+b
        return f

    def sq_fit(x,a,b,c):
        '''
        quadratic function (y=ax**2+b*x+c)
        '''
        f = (a*(x**2))+(b*x)+c
        return f

    def subhalo_url_constructor(snapID, subID):
        url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}".format(str(snapID),str(subID))
        return url


class subhalo_tree:
    def __init__(self, snapID, subID):
        r'''
        Load known snapshot and subhalo ID for trees tracing request 
        '''
        
        self.snapID = snapID #initial known value 
        self.subID = subID  # initial known value 
        
        self.primesuburl = UTILITY.subhalo_url_constructor(self.snapID, self.subID)
        self.primesubhalo = UTILITY.get(self.primesuburl)
        sub = self.primesubhalo
        #'''
        #self.mpb1 = UTILITY.treeget(sub['trees']['sublink_mpb'])
        self.mpb2 = UTILITY.treeget(sub['trees']['lhalotree_mpb'])
        '''
        try:
            self.mpb1 = 'trees/sublink_mpb_{}.hdf5'.format(self.subID)
            self.mpb2 = 'trees/lhalotree_mpb_{}.hdf5'.format(self.subID)
        except FileNotFoundError:
            self.mpb1 = UTILITY.treeget(sub['trees']['sublink_mpb'])
            self.mpb2 = UTILITY.treeget(sub['trees']['lhalotree_mpb'])
        except IOError:
            self.mpb1 = UTILITY.treeget(sub['trees']['sublink_mpb'])
            self.mpb2 = UTILITY.treeget(sub['trees']['lhalotree_mpb'])
        '''
        
    def idtrace(self):
        with h5py.File(self.mpb2,'r') as f:
            snapnums = f['SnapNum'][:]
            subid = f['SubhaloNumber'][:]
        snapnum = list(snapnums); subid = list(subid)
        snapnum.reverse();subid.reverse()
        df_id = pd.DataFrame({
            "snapshot": snapnum,
            "id": subid
        })
        self.snapnum = snapnum
        self.subids = subid
        self.id_frame = df_id.copy()
        #print (self.snapnum)
        #print(len(self.snapnum))        
        return df_id
    
    def subhalodata(self,i):
        url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/".format(self.snapnum[i], self.subids[i])
        sub = UTILITY.get(url)
        met = sub['gasmetallicity']+sub['starmetallicity']
        mass = sub['mass_log_msun']
        sfr = sub['sfr']
        print("snap{} done!".format(self.snapnum[i]))
        return (met,mass,sfr,self.subids[i],self.snapnum[i])
    
    def MSDATAGET(self):
        returns = Parallel(n_jobs= 2)(delayed(self.subhalodata)(i) for i in range(len(self.snapnum)))
        df=pd.DataFrame(returns,columns=['met','mass','sfr','id','snapshot'])
        self.maindf = df
        return df
    
    def MSPLOT(self):
        snapshots = [1,10,21,33,40,50,67,78,91,99]
        linesnaps = [1,10,21,33,40,50,67,78,91]
        linedf = self.maindf[self.maindf['snapshot'].isin(linesnaps)]
        self.maindf = self.maindf[self.maindf['snapshot'].isin(snapshots)]
        print(self.maindf)
        basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'
        fields = ['SubhaloMass','SubhaloSFRinRad','SubhaloSFR','SubhaloGasMetallicitySfrWeighted']
        subhalos99 = il.groupcat.loadSubhalos(basePath,67,fields=fields)
        mass =subhalos99['SubhaloMass'] * 1e10 / 0.704
        sfr = list((subhalos99['SubhaloSFR']))
        plt.figure(figsize=(20,12))
        plt.yscale('log')
        plt.title("evolution of largest progenitor to subhalo 21 at z=0: ")
        plt.plot(np.log10(mass),sfr, 'g+',label = 'Snapshot 99 subhalos', zorder=1)
        plt.plot(linedf['mass'],linedf['sfr'], 'r--',label = 'path of progenitors to subhalo 21',zorder=2)
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
        plt.savefig('msev_{}.png'.format(self.subID))
        plt.close()
        
    def plot3d(self):
        baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/"
        search_q = "?limit=17553&sfr__gt=0.0"
        url = baseurl+ search_q
        subs = UTILITY.get(url)
        ids = [subs['results'][i]['id'] for i in range(subs['count'])]
        mass = []
        sfr = []
        urls = []
        df = self.maindf
        for i,id in enumerate(ids):
            mass.append(subs['results'][i]['mass_log_msun'])
            sfr.append(subs['results'][i]['sfr'])
            urls.append(subs['results'][i]['url'])
        plt.style.use('_mpl-gallery')
        
        fig = plt.figure()
        ax = plt.axes(projection="3d")  
        ax.plot3D(mass,np.log10(sfr),zs=0, zdir='z', c='g', marker = '+', linestyle='None')
        ax.scatter3D(df['mass'],np.log10(df['sfr']),df['snapshot'],c=df['id'],cmap = 'magma')
        ax.plot3D(df['mass'],np.log10(df['sfr']),df['snapshot'],c='r',linestyle = 'dashed')
        fig.savefig('3d19.png')
        plt.show()
        
#
#Highest Resolution Snapshots 
#Corresponds to z=[4.01, 2.00, 1.00, 0.50, 0.30, 0.10, 0.00]
#



def subhalo_traces(i):
    try:
        test = subhalo_tree(99,i)
        test.idtrace()
        test.MSDATAGET()
        test.MSPLOT()
        return print("done for subhalo{}".format(i))
    except OSError as e:
        return print(e)
    '''
    except TypeError as e:
        return print(e)
    except IndexError as e:
        return print(e)
    except ValueError as e:
        return print(e)
    '''
subhalo_traces(106)
#100,101,108,115,120,121,122
#returns = Parallel(n_jobs = 20)(delayed(subhalo_traces)(i) for i in ids)