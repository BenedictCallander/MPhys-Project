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

class cutsub:
    def __init__(self,subID,snapID,simID,primeID):
        
        self.subID = subID
        self.snapID = snapID
        self.simID = simID
        self.primeID = primeID
        redurl = 'http://www.tng-project.org/api/'+str(simID)+'/snapshots/'+str(snapID)
        utility = UTILITY.get(redurl)
        redshift = utility['redshift']
        self.redshift =redshift
        hubble = 0.7
        scalefac = 1./(1.+redshift) #calculate scale factor

        #
        #investigate particle level subhalo data 
        #
        # begin by downloading subhalo cutout       
        self.subURL = 'http://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/'.format(str(snapID),str(subID))
        subhalo = UTILITY.get(self.subURL)
        #
        #obtain global properties
        #
        self.mass = subhalo['mass_log_msun']
        self.tot_sfr = subhalo['sfr']
        self.tot_met = subhalo['gasmetallicity']
        self.Rhalf = subhalo['halfmassrad']
        self.crit_dist = 5* self.Rhalf
        self.stellarphotometrics = subhalo['stellarphotometricsrad']
        
        
        self.centre = np.array([subhalo['pos_x'],subhalo['pos_y'],subhalo['pos_z']])
        
        cutout ="files/historycutouts/evdir_{}/cutout_{}.hdf5".format(self.primeID,self.subID)
        with h5py.File(cutout,'r') as f:
            sfr = f['PartType0']['StarFormationRate'][:]
            co_ords = f['PartType0']['Coordinates'][:]
            hcoldgas  = np.where( (sfr > 0.0))[0]
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
        f=df.sample(frac=0.1,replace=False)

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
    
    def linearfit(self, dfin):
        dfin.sort_values(by='rad',inplace=True)
        popt,pcov = curve_fit(UTILITY.linear_fit, dfin['rad'],np.log10(dfin['met'])+12,sigma=1/dfin['sfr'],absolute_sigma=True)
        med_data = medfilt(np.log10(dfin['met'])+12,kernel_size = 21)

        plt.figure(figsize=(20,12))
        plt.title("Metgrad for {} - snap{} (linked to sub{}_snap99)".format(self.subID,self.snapID,subhaloid))
        plt.plot(dfin['rad'], med_data, 'r-')
        plt.plot(dfin['rad'], UTILITY.linear_fit(dfin['rad'],*popt))
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}O/H)$ (SFR Normalised)")
        plt.ylim(8,11)
        filename = 'historypng/snap{}_progenitorto_{}png'.format(self.snapID,subhaloid)
        plt.savefig(filename)
        plt.close()
        #print("gradient {}".format(popt[0]))
        return popt[0]
    
    def piecewise(self,dfin,breakpoint):
        df = dfin.copy()
        df.sort_values(by="rad",inplace = True)
        x0 = np.array([min(df['rad']), breakpoint, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met']),weights=1/df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        slope1 = my_pwlf.slopes[0]
        slope2 = my_pwlf.slopes[1]
        #print("slopes are inner: {} and outer:{}".format(slope1,slope2))
        '''
        med_data1 = medfilt((12+np.log10(df['met'])), kernel_size=11)
        xHat = np.linspace(min(df['rad']), max(df['rad']), num=10000)
        yHat = my_pwlf.predict(xHat)
        plt.figure(figsize=(20,12))
        plt.plot(df['rad'], med_data1, 'b--')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        #filename = 'histbrfit/single/{}_snap_sub={}.png'.format(self.snapID, self.subID)
        plt.savefig("testing{}.png".format(self.subID))
        plt.close()
        '''
        return (slope1,slope2)
    
    def doublepiecewise(self,dfin,breakpoint1,breakpoint2):
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
        plt.figure(figsize=(20,12))
        plt.plot(df['rad'], med_data1, 'b--')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        filename = 'histbrfit/double/{}_sub_{}_doublebreak.png'.format(self.snapID, self.subID, self.snapID)
        plt.savefig(filename)
        plt.close()

df = pd.read_csv("all.csv")
subs = list(df['subhalo'])
snaps = list(df['snapshot'])
prime = list(df['primesub'])
l1 = []; l2= []
def dofunc(i):
    try:
        sub = cutsub(subs[i],snaps[i],'TNG50-1',prime[i])
        sub.align_dfgen()
        df = sub.filter()
        slope1,slope2 = sub.piecewise(df,3)
        subID = subs[i];snapID = snaps[i]; primeID = prime[i]
        print("done for sub {} at snap {}".format(subID,snapID))
        return (slope1,slope2,subID,snapID,primeID)
    except OSError as e:
        return print(e)
    except TypeError as e:
        return print(e)
    except IndexError as e:
        return print(e)
    except ValueError as e:
        return print(e)

returns = Parallel(n_jobs=20)(delayed(dofunc)(i) for i in range(len(subs)))
dfout = pd.DataFrame(returns, columns = ['slope1', 'slope2', 'subhalo', 'snapshot', 'primesub'])
dfout.to_csv("all2.csv")