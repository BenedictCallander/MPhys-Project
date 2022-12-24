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
    
    def median_filter(self,dfin):
        df = dfin.copy()
        #print(df)
        df.sort_values(by='rad',inplace=True)
        ranges = np.arange(df.rad.min() - 0.1, df.rad.max() + 0.1, 0.1)  
        groups = df.groupby(pd.cut(df.rad, ranges))
        medianout = groups.median()
        medianout = medianout.dropna()
        return medianout
    
    def linearfit(self, dfin):
        dfin.sort_values(by='rad',inplace=True)
        popt,pcov = curve_fit(UTILITY.linear_fit, dfin['rad'],np.log10(dfin['met'])+12,sigma=1/dfin['sfr'],absolute_sigma=True)
        return popt[0]
    
    def piecewise(self,dfin,breakpoint):
        df = dfin.copy()
        x0 = np.array([df.rad.min(), breakpoint, df.rad.max()])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met']),weights=df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        slope1 = my_pwlf.slopes[0]
        slope2 = my_pwlf.slopes[1]
        return (slope1,slope2)
    
    def getbreaks(self,dfin):
        df = dfin.copy()
        x = list(df['rad'])
        y = list(df['met'])
        piecewise_model = pwlf.PiecewiseLinFit(x, y)
        num_breaks = 1
        piecewise_model.fit(num_breaks)
        breaks = piecewise_model.fit_guess([6.0])
        return breaks[1]
    


class dodirectory:
    def __init__(self,primeID):
        self.primeID = primeID
    def getlist(self):
        listpath = "files/historycutouts/evdir_{}/treedata_{}.csv".format(self.primeID,self.primeID)
        df = pd.read_csv(listpath)
        snapshots = list(df['snapshots'])
        subhalos = list(df['subhalos']) 
        return snapshots,subhalos

def dosingle(sub,snap,prime):
    try:
        subhalo = cutsub(sub,snap,'TNG50-1', prime)
        subhalo.align_dfgen()
        df= subhalo.filter()
        df = subhalo.median_filter(df)
        subID = sub
        snapID = snap
        breakrad = subhalo.getbreaks(df)
        print("done for subhalo {} snapshot {}".format(sub,snap))
        return (subID, snapID,breakrad)
    except OSError as e:
        return print(e)
    except TypeError as e:
        return print(e)
    except KeyError as e:
        return print(e)
    except IndexError as e:
        return print(e)
    except ValueError as e:
        return print(e)
    
def dodir(i):
    try:
        data = dodirectory(i)
        snapshots,subhalos = data.getlist()
        subs = [];snaps=[];breakrads = []
        for j in range(len(snapshots)):
            subID,snapID,breakrad = dosingle(subhalos[j],snapshots[j],i)
            subs.append(subID);snaps.append(snapID)
            breakrads.append(breakrad)
        df = pd.DataFrame({
            'subhalo':subs,
            'snapshot':snaps,
            'breakrad':breakrad
        })
        #returns = Parallel(n_jobs=4)(delayed(dosingle)(subhalos[j],snapshots[j],i)for j in range(4))
        #df = pd.DataFrame(returns,columns = ['subhalo','snapshot','slope1','slope2'])
        fpath = "files/historycutouts/evdir_{}/slopebreakpoint{}.csv".format(i,i)
        df.to_csv(fpath)
        return print("done for descendant {}".format(i))
    except OSError as e:
        return print(e)
    except TypeError as e:
        return print(e)
    except KeyError as e:
        return print(e)
    except IndexError as e:
        return print(e)
    except ValueError as e:
        return print(e)
    
dfyay = pd.read_csv("csv/traceids.csv")
ids = list(dfyay['id'])
#for i in ids:
#    dodir(i)
returns = Parallel(n_jobs=25)(delayed(dodir)(i) for i in ids)

#208568
'''
@numba.jit
def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def fit_piecewise_linear(xdata, ydata):
    p, e = curve_fit(piecewise_linear, xdata, ydata, (5, 8, 1, -1))
    x0, y0, k1, k2 = p
    return k1, k2

xdata = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
ydata = [1, 4, 6, 7, 8, 9, 9, 8, 6, 4]

k1, k2 = fit_piecewise_linear(xdata, ydata)
print(k1, k2)
'''