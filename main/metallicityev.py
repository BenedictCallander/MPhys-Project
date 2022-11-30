# Metallicityev.py 
# \-> script containing Classes subsequent functions to study the Metallicity evolution of a subhalo's metallicity gradient through the IllustrisTNG snapshots 
# Created:17/11/2022 
# Author: Benedict Callander 
# Respository https://github.com/btcallander/MPhys-Project (private)
#

#Plotting, numerical functions and dataset manipulation
import matplotlib.pyplot as plt
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
import BCUTILS

#runtime calculation 
import time

#Computational functions - simultaneous calculations to make use of multi-core CPU
from joblib import Parallel, delayed

#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter

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
    

class subhaloev:
    def __init__(self, startID, startSnap, simID,):
        self.startID = startID
        self.startSN = startSnap
        self.simID = simID
        
        self.start_url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}".format(str(self.startSN), str(self.startID))

    def fetchtree(self):
        subhalo = UTILITY.get(self.start_url)
        self.mpb2 = UTILITY.get(subhalo['trees']['lhalotree_mpb'])
        
    def fetchIDS(self):
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
        
    def get_history(self,i):
        url = "https://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/".format(self.snapnum[i], self.subids[i])
        sub = UTILITY.get(url)
        met = sub['gasmetallicity']+sub['starmetallicity']
        mass = sub['mass_log_msun']
        sfr = sub['sfr']
        print("snap{} done!".format(self.snapnum[i]))
        return (met,mass,sfr,self.subids[i],self.snapnum[i])
    
    def history_filegen(self, filepath):
        returns = Parallel(n_jobs= 4)(delayed(self.get_history)(i) for i in range(1,99))
        df=pd.DataFrame(returns,columns=['met','mass','sfr','id','snapshot'])
        df.to_csv(filepath)
                    
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
        
class subhalo_cut:
    def __init__(self, subID, snapID, simID):
        self.subID = subID
        self.snapID = snapID
        self.simID = simID
        
        self.baseurl = "https://www.tng-project.org/api/{}/snapshots/{}/subhalos/{}/".format(str(self.simID),self.snapID, self.subID)
        
    def get_cutout(self):
        baseurl = self.baseurl
        subhalo = UTILITY.get(self.baseurl)
        cutout_request = {'gas':'Coordinates,Masses,GFM_Metallicity,StarFormationRate'}
        cutout = UTILITY.get(baseurl+"cutout.hdf5",cutout_request)
        
        with h5py.File(cutout,'r') as f:
            x = f['PartType0']['Coordinates'][:,0] - subhalo['pos_x']
            y = f['PartType0']['Coordinates'][:,1] - subhalo['pos_y']
            z = f['PartType0']['Coordinates'][:,2] - subhalo['pos_z']
            sfr = f['PartType0']['StarformationRate']
            metallicity = f['PartType0']['GFM_Metallicity']
        df = pd.DataFrame({
            "x":x,
            "y":y,
            "z":abs(z),
            "sfr":sfr,
            "metallicity":12+np.log10(metallicity)
        })
        
        zmax = 0.1*subhalo['stellarphotometricsrad']
        dfnew = df[df['z']<zmax]
        rad = np.sqrt(df['x']**2+df['y']**2)
        dfnew.insert(4,"rad", rad, True)
        dfnew.rad = 10*((dfnew.rad-dfnew.rad.min())/(dfnew.rad.max()-dfnew.rad.min()))
        
        self.df = dfnew
        
    def met_gradient(self):
        df = self.df
        
        medfitt = medfilt(df['met'], kernel_size=3)
        popt,pcov = curve_fit(UTILITY.linear_fit,df['rad'],df['metallicity'], sigma = 1/df['sfr'], absolute_sigma=True)
        raddata = np.linspace(0,10,100)
        plt.figure(figsize=(20,12))
        plt.plot(df['rad'], medfitt,'b-')
        plt.plot(df['rad'], UTILITY.linear_fit(df['rad'],*popt),'r--')
        plt.xlabel("radial distance")
        plt.ylabel("12+$log_{10}$(O/H)")
        filename = "historypng/metgrad_snap_{}_sub_{}.png".format(self.snapID,self.subID)
        plt.savefig(filename)
        plt.close()
        
        
subhalo = subhaloev(11,99,'TNG50-1')
        
        
        
        
        
        
            
        