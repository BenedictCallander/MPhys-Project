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

#treepath = "files/historycutouts/evdir_{}/treedata_{}.csv".format(start,start)

class cutout:
    def __init__(self, start, subID,snapID,rhalf,stellarphotometrics):
        self.startsub = start
        self.subID = subID
        self.snapID = snapID
        #
        #
        #
        self.Rhalf = rhalf
        self.stellarphotometrics = stellarphotometrics     
        cutpath = "files/historycutouts/evdir_{}/cutout_{}.hdf5".format(start,subID)
        with h5py.File(cutpath,'r') as f:
            sfr = f['PartType0']['StarFormationRate'][:]
            co_ords = f['PartType0']['Coordinates'][:]
            hcoldgas  = np.where( (sfr > 0.0))[0]
            self.pgas_coo = f['PartType0']['Coordinates'][hcoldgas]
            self.pgas_m = f['PartType0']['Masses'][hcoldgas]
            self.pgas_vel = f['PartType0']['Velocities'][hcoldgas]
            self.pgas_met = f['PartType0']['GFM_Metallicity'][hcoldgas]
            self.pgas_sfr = f['PartType0']['StarFormationRate'][hcoldgas]
        #
        self.test = len(hcoldgas)
        print(self.test)
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
        self.df=df[df['rad']<5*self.Rhalf]
        
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
        x0 = np.array([min(df['rad']), breakpoint1,breakpoint2, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met']),weights=1/df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        slope1 = my_pwlf.slopes[0]
        slope2 = my_pwlf.slopes[1]
        slope3 = my_pwlf.slopes[2]
        
        print("slopes are inner: {} middle:{} and outer:{}".format(slope1,slope2,slope3))
        '''
        xHat = np.linspace(min(df['rad']), max(df['rad']), num=10000)
        yHat = my_pwlf.predict(xHat)
        plt.figure(figsize=(20,12))
        plt.plot(df['rad'], med_data1, 'b--')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        filename = 'files/images/metallicityhistory/{}_sub_{}.png'.format(self.primesub, self.snapID, self.subID)
        plt.savefig(filename)
        plt.close()
        '''
        

class dirrun:
    def __init__(self,descendant):
        self.primesub = descendant
    
    def getlist(self):
        '''
        f['list1'].attrs['keyword'] = 'mass'
            f['list2'].attrs['keyword'] = 'sfr'
            f['list3'].attrs['keyword'] = 'met'
            f['list4'].attrs['keyword'] = 'Rhalf'
            f['list5'].attrs['keyword'] = 'stellarph'
            f['list6'].attrs['keyword'] = 'snaps'
            f['list7'].attrs['keyword'] = 'subs'
        '''
        datapath = "files/historycutouts/evdir_{}/subdata_{}.hdf5".format(self.primesub,self.primesub)
        with h5py.File(datapath,'r') as f:
            snapshots = list(f['list6'][:])
            subhalos = list(f['list7'][:])
            Rhalf = list(f['list4'][:])
            stellarphotometrics = list(f['list5'][:])
            
            return snapshots,subhalos,Rhalf,stellarphotometrics
    
    def dir(self):
        datapath = "files/historycutouts/evdir_{}/subdata_{}.hdf5".format(self.primesub,self.primesub)
        with h5py.File(datapath,'r') as f:
            snapshots = list(f['list6'][:])
            subhalos = list(f['list7'][:])
            Rhalf = list(f['list4'][:])
            stellarphotometrics = list(f['list5'][:])
            l1 = []; l2=[];l3 = []
            for i in range(len(snapshots)):
                sub = cutout(self.primesub, subhalos[i],snapshots[i], Rhalf[i], stellarphotometrics[i])
                sub.align_dfgen()
                df = sub.filter()
                slope1,slope2,slope3 = sub.doublepiecewise(df,3,8)
                l1.append(slope1); l2.append(slope2);l3.append(slope3)
            print(l1,l2,l3)
            
            
test = dirrun(63870)
test.dir()

end = time.time()
print("runtime = {}".format(end-start))