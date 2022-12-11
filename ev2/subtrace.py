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

#hdf5 binary file manipulation
import h5py

#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il


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


class subhistory:
    def __init__(self,startID, startSN):
        self.startID = startID
        self.startsnap = startSN
        
        self.start_url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}".format(str(startSN),str(startID))
        self.start_sub = UTILITY.get(self.start_url)
        self.check = self.start_sub['related']['sublink_progenitor']
        #self.mpb = UTILITY.get(self.start_sub['trees']['sublink_mpb'])
        #self.mdb = UTILITY.get(self.start_sub['trees']['sublink_mdb'])
        
        self.mpb = "files/trees/sublink_mpb_{}.hdf5".format(startID)
        #self.mdb = "trees1/desc/sublink_mdb_{}.hdf5".format(startID)
        
        # collect all IDS 
        with h5py.File(self.mpb,'r') as f:
            presnap = list(f['SnapNum'][:])
            preID= list (f['SubfindID'][:])
        '''  
        with h5py.File(self.mdb,'r') as x:
            postsnap = list(x['SnapNum'][:])
            postID = list(x['SubfindID'][:])
        postsnap.reverse() ; postID.reverse()
        presnap.extend(postsnap)
        preID.extend(postID)
        
        '''
        presnap.reverse(); preID.reverse()
              
        self.snapnum = presnap
        self.subids = preID
        
    def subhalodata(self,i):
        url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/".format(self.snapnum[i], self.subids[i])
        sub = UTILITY.get(url)
        met = sub['gasmetallicity']+sub['starmetallicity']
        mass = sub['mass_log_msun']
        sfr = sub['sfr']
        print("snap{} done!".format(self.snapnum[i]))
        return (met,mass,sfr,self.subids[i],self.snapnum[i])
    
    def MSDATAGET(self):
        returns = Parallel(n_jobs= 20)(delayed(self.subhalodata)(i) for i in range(len(self.snapnum)))
        df=pd.DataFrame(returns,columns=['met','mass','sfr','id','snapshot'])
        self.maindf = df
        return df
    
    def MSPLOT(self):
        snapshots = [1,5,6,10,21,33,40,50,52,55,67,78,91,99]
        linesnaps = [1,5,6,10,21,33,40,50,52,55,67,78,91]
        linedf = self.maindf[self.maindf['snapshot'].isin(linesnaps)]
        self.maindf = self.maindf[self.maindf['snapshot'].isin(snapshots)]
        print(self.maindf)
        snap = 67
        basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'
        fields = ['SubhaloMass','SubhaloSFRinRad','SubhaloSFR','SubhaloGasMetallicitySfrWeighted']
        subhalos99 = il.groupcat.loadSubhalos(basePath,snap,fields=fields)
        mass =subhalos99['SubhaloMass'] * 1e10 / 0.704
        sfr = list((subhalos99['SubhaloSFR']))
        plt.figure(figsize=(20,12))
        plt.yscale('log')
        plt.title("evolution of largest progenitor to subhalo {} at z=0: ".format(self.startID))
        plt.plot(np.log10(mass),sfr, 'g+',label = 'Snapshot 99 subhalos', zorder=1)
        plt.plot(linedf['mass'],linedf['sfr'], 'r--',label = 'path of progenitors to subhalo {}'.format(self.startID),zorder=2)
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
        plt.savefig('msev_{}_bg{}.png'.format(self.startID,snap))
        plt.close()

def runscript(i):
    sub = subhistory(i,99)
    sub.MSDATAGET()
    sub.MSPLOT()
    return print("subhalo {} done".format(i))

runscript(37)
#df= pd.read_csv('traceids.csv')
#ids = list(df['id'])

#returns = Parallel(n_jobs= 2)(delayed(runscript)(i) for i in ids)


'''
for i in ids:
    runscript(i)
    print("subhalo {} done".format(i))

IDS
[199, 226, 250, 261, 265, 282, 287, 288, 302, 343, 351, 54120, 54123, 92900, 92906, 92943, 92954, 105165, 115190, 115194, 160834, 160854,
166626, 166632, 166641, 173159, 173161, 173163, 173168, 173170, 173174, 173175, 173183, 173185, 173190, 173199, 190064, 237953, 237957,
237961, 237966, 237967, 237969, 242568, 244821, 244822, 244827, 244829, 244835, 261422, 261426, 261427, 261432, 261435, 261438, 261439, 262794,
266542, 268475, 268476, 268477, 268481, 279688, 281000, 281001, 281002, 281003, 281008, 282221, 282225, 282226, 289162, 289165, 289167, 289171,
296873, 308513, 308515, 308523, 317887, 317888, 317889, 324131, 324135, 327105, 331162, 331165, 334220, 334221, 334222, 334223, 334224, 334225,
334226, 340064, 348016, 353213, 354844, 358981, 359436, 359437, 359438, 359439, 374420, 374421, 374422, 377803, 379570, 382231, 385419, 385420, 388001,
388002, 388580, 388582, 388583, 404820, 438731, 442224]

'''