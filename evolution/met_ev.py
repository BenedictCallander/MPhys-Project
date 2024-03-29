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


class cutout_subhalo:
    def __init__(self,snapshot,subhalo,descendant):
        filedir = "files/binary/historycutouts/evdir_{}/".format(descendant)
        self.subID = subhalo
        self.snapID = snapshot
        self.startsub = descendant
        snapshots = [21,33,50,67,78,91,99]
        redshifts= [4.01, 2,1,0.5,0.3,0.1,0.0]
        index = snapshots.index(self.snapID)
        redshift= redshifts[index]
        self.redshift = redshift    
        scalefac = 1./(1.+redshift) #calculate scale factor
  
        #
        #Query API for global data 
        #
        self.subURL = 'http://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/'.format(str(self.snapID),str(self.subID))
        subhalo = UTILITY.get(self.subURL)
        
        self.mass = subhalo['mass_log_msun']
        self.tot_sfr = subhalo['sfr']
        self.tot_met = subhalo['gasmetallicity']
        self.Rhalf = subhalo['halfmassrad']
        self.crit_dist = 5* self.Rhalf
        self.stellarphotometrics = subhalo['stellarphotometricsrad']
        self.centre = np.array([subhalo['pos_x'],subhalo['pos_y'],subhalo['pos_z']])

        cutout = filedir+"cutout_{}.hdf5".format(self.subID)
        with h5py.File(cutout,'r') as f:
            sfr = f['PartType0']['StarFormationRate'][:]
            co_ords = f['PartType0']['Coordinates'][:]
            hcoldgas  = np.where( (sfr > 0.0))[0]
            #print(sfr)
            #print(hcoldgas)
            self.pgas_coo = f['PartType0']['Coordinates'][hcoldgas]
            self.pgas_m = f['PartType0']['Masses'][hcoldgas]
            self.pgas_vel = f['PartType0']['Velocities'][hcoldgas]
            self.pgas_met = f['PartType0']['GFM_Metallicity'][hcoldgas]
            self.pgas_sfr = f['PartType0']['StarFormationRate'][hcoldgas]
        self.test = len(hcoldgas)

        self.pgas_coo -= self.centre[None,:]
        
    def align_dfgen(self):
        _coo = np.copy(self.pgas_coo)
        _vel = np.copy(self.pgas_vel)
        _m = np.copy(self.pgas_m)
        
        self.ang_mom_3D = np.sum(_m[:,None,]*np.cross(_coo,_vel),axis=0)
        self.ang_mom = self.ang_mom_3D/ np.sum(_m)

        j=self.ang_mom/np.linalg.norm(self.ang_mom)
        #normalised specific angular momentum 
        
        x = np.array([1,2,3])
        x = x-(x.dot(j)*j) #make x orthogonal to j
        
        x/= np.linalg.norm(x) # normalise
        
        y = np.cross(j,x)#create 3rd vector - orth to x,j
        
        
        A = (x,y,j) # transformation matrix
        
        self.pgas_coo=np.dot(A,self.pgas_coo.T).T # change co-ordinates
        self.pgas_vel = np.dot(A,self.pgas_vel.T).T
        self.radial = np.sqrt((self.pgas_coo[:,0]**2)+(self.pgas_coo[:,1]**2))
        
        df = pd.DataFrame({
            "x":self.pgas_coo[:,0],
            "y":self.pgas_coo[:,1],
            "z":self.pgas_coo[:,2],
            "rad": self.radial,
            "mass":self.pgas_m,
            "met":(self.pgas_met),
            "sfr":self.pgas_sfr
            })
        df=df[df['rad']<5*self.Rhalf]
        #print(df)
        self.df = df
        
    def filter(self):
        df = self.df.copy()
        spr = self.stellarphotometrics
        z_max = 0.1*spr
        df = df[df['z']<z_max]
        df.rad =10*((df.rad-df.rad.min())/(df.rad.max()-df.rad.min()))
        
        self.df = df
        return df
    
    def doublepiecewise(self,dfin,breakpoint1,breakpoint2):
        dfin.sort_values(by='rad',inplace = True)
        df = dfin.copy()
        df.sort_values(by="rad",inplace = True)
        med_data1 = medfilt((12+np.log10(df['met'])), kernel_size=11)
        x0 = np.array([min(df['rad']), breakpoint1,breakpoint2, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met']),weights=1/df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        slope1 = my_pwlf.slopes[0]
        slope2 = my_pwlf.slopes[1]
        slope3 = my_pwlf.slopes[2]
        
        print("slopes are inner: {} middle:{} and outer:{}".format(slope1,slope2,slope3))
        xHat = np.linspace(min(df['rad']), max(df['rad']), num=10000)
        yHat = my_pwlf.predict(xHat)
        '''
        plt.figure(figsize=(20,12))
        plt.plot(df['rad'], med_data1, 'b--')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        filename = 'files/images/metallicityhistory/99progenitors/progenitorsto_{}/{}_sub_{}_doublebreak.png'.format(self.startsub, self.snapID, self.subID)
        plt.savefig(filename)
        plt.close()
        '''
        
        return (slope1,slope2,slope3)



class metallicity_evolution:
    def __init__(self,descendant):
        self.startsub = descendant
        
        filedir = "files/binary/historycutouts/evdir_{}/".format(self.startsub)
        df = pd.read_csv(filedir+"treedata_{}.csv".format(self.startsub))
        target_snaps = list(df['snapshots'])
        target_subs = list(df['subhalos'])
        
        self.target_snaps = target_snaps
        self.target_subs = target_subs
    
    def individual_subhalo(self,i):
        snapshot = self.target_snaps[i]
        subhalo = self.target_subs[i]
        subhalo = cutout_subhalo(snapshot, subhalo, self.startsub)
        subhalo.align_dfgen()
        dfg = subhalo.filter()
        slope1,slope2,slope3= subhalo.doublepiecewise(dfg,3,8)
        self.slope1 = slope1
        self.slope2 = slope2
        self.slope3 = slope3
        return (slope1,slope2,slope3)
    
    def history_trace(self):
        fills = [self.startsub,self.startsub,self.startsub,self.startsub,self.startsub,self.startsub]
        snapshots = self.target_snaps
        subhalos = self.target_snaps
        
        inners = []; middles =  [] ; outers= []
        for i in range(6):
            slope1,slope2,slope3 = metallicity_evolution.individual_subhalo(self,i)
            inners.append(slope1)
            middles.append(slope2)
            outers.append(slope3)
        df = pd.DataFrame({
            "progenitor": fills,
            "snapshot":snapshots,
            "subhalo": subhalos,
            "slope1": inners,
            "slope2": middles,
            "slope3": outers
        })
        fpath = "files/binary/historycutouts/evdir_{}/slopehistory_{}".format(self.startsub,self.startsub)
        df.to_csv(fpath)
             

df = pd.read_csv("traceids.csv")
ids = list(df['id'])

def main(i):
    try:
        history = metallicity_evolution(i)
        history.history_trace()
        return print("done for {}".format(i))
    except OSError as e:
        return print(e)
    except TypeError as e:
        return print(e)
    except IndexError as e:
        return print(e)
    except ValueError as e:
        return print(e)
    
returns = Parallel(n_jobs= 20)(delayed(main)(i) for i in ids)
