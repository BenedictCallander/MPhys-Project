from itertools import groupby
import logging
from random import random
from re import sub # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il
from scipy.signal import medfilt
from scipy.optimize import curve_fit
from joblib import Parallel, delayed
import os
from scipy.signal import savgol_filter
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
#basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output'
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

def linear_fit(a,x,b):
    f = (a*x)+b
    return f
def sq_fit(x,a,b,c):
    f = (a*(x**2))+(b*x)+c
    return f


class galaxy:
    def __init__(self,simID,snapID,subID,SFRYN):
        np.seterr(divide='ignore', invalid='ignore')
        #object creation requires only 3 input parameters for galaxy selection 
        self.simID = simID #simulation used (TNG100/TNG50)
        self.snapID = snapID #snapshot of study - can be either number of z= (z= required in parenthesis )
        self.subID = subID #subhalo ID - could be iterated 
        #basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output' #path to simulation data on Noether
        basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'
        baseurl = 'http://www.tng-project.org/api/'+str(simID)+'/snapshots/'+str(snapID) # API url for simulation/snapshot information
        
        hubble =0.7
        
        snap_util =get(baseurl) #use api to easily read snapshot information
        redshift = 0
        self.redshift =redshift  # store redshift value as attribute
        
        scalefac = 1./(1.+redshift) #calculate scale factor
        
        #
        # Read Subhalo level info 
        #

        ptNumGas = il.snapshot.partTypeNum('gas') #determine index designation for each particle type
        ptNumStars = il.snapshot.partTypeNum('stars')
        #pull all data for specific subhalo 
        all_fields= il.groupcat.loadSingle(basePath, snapID, subhaloID = subID)
        self.test=all_fields['SubhaloMassInRadType'][ptNumGas]
        
        self.lMgas  = np.log10( all_fields['SubhaloMassInRadType'][ptNumGas]/hubble ) + 10.
        self.lMstar = np.log10( all_fields['SubhaloMassInRadType'][ptNumStars]/hubble ) + 10.
        # Coordinate of particle with minimum binding energy (converted from ckpc/h to kpc)
        self.centre = all_fields['SubhaloPos']/hubble / (1. + redshift)  # 3-element array [units: proper kpc]
        # Adopt the 3D half-stellar-mass radius
        self.Rhalf  = all_fields['SubhaloHalfmassRadType'][ptNumStars]/hubble / (1. + redshift)  
        self.stellarphotometricsrad = all_fields['SubhaloStellarPhotometricsRad']
        # [units: proper kpc] (quantified in 3D)
        
        # Load all the relevant particle level info
        gas = il.snapshot.loadSubhalo(basePath, snapID, subID, 'gas', fields=['Coordinates', 'Masses','Density','Velocities', 'StarFormationRate','GFM_Metallicity'])
        # dimensions and units (see https://www.tng-project.org/data/docs/specifications/#parttype0):
        # Coordinates (N,3) ckpc/h   where ckps stands for co-moving kpc
        # Masses      (N)   10**10 Msun/h
        # Velocities  (N,3) km sqrt(scalefac)        # We convert these to pkpc (proper kpc), Msun and km/s, respectively
        crit_dist = 5 * self.Rhalf #30. # proper kpc
        self.crit_dist = crit_dist
        if (SFRYN=='Y'):
            hcoldgas  = np.where( (gas['StarFormationRate'] > 0.) & (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        elif(SFRYN=='N'):
            hcoldgas  = (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2)
        self.pgas_coo   = gas['Coordinates'][hcoldgas]/hubble / (1. + redshift)
        self.pgas_m     = gas['Masses'][hcoldgas] * 10**10 / hubble
        self.pgas_vel   = (gas['Velocities'][hcoldgas] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.conv_kms2kpcyr = (3.1558 / 3.08568) * 10**(-9)
        self.pgas_vel   = self.pgas_vel * self.conv_kms2kpcyr    #Convert to kpc/yr
        self.pgas_sfr   = gas['StarFormationRate'][hcoldgas]
        self.pgas_met   =gas['GFM_Metallicity'][hcoldgas]
        self.pgas_dens = gas['Density'][hcoldgas]
        #print(all_fields.keys())
        # Load all stellar particle data
        stars = il.snapshot.loadSubhalo(basePath, snapID, subID, 'stars', fields=['Coordinates', 'Masses', 'Velocities','GFM_Metallicity' ])
        hstar = np.where( (np.sum((stars['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        self.pstar_coo   = stars['Coordinates'][hstar]/hubble / (1. + redshift)
        self.pstar_m     = stars['Masses'][hstar] * 10**10 / hubble
        self.pstar_vel   = (stars['Velocities'][hstar] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.pstar_vel   = self.pstar_vel * self.conv_kms2kpcyr
        self.pstar_met = stars['GFM_Metallicity'][hstar]


    def galcen(self):
        self.pgas_coo -= self.centre[None,:]
        self.pstar_coo -= self.centre[None,:]
        
    def ang_mom_align(self, type):
        if (type=='gas'):
            _coo = np.copy(self.pgas_coo)
            _vel = np.copy(self.pgas_vel)
            _m = np.copy(self.pgas_m)
        elif(type=='stars'):
            _coo =np.copy(self.pstar_coo)
            _vel = np.copy(self.pstar_vel)
            _m = np.copy(self.pstar_m)
        # calc angular momentum based on particle type 
        
        self.ang_mom_3D = np.sum(_m[:,None,]*np.cross(_coo,_vel), axis = 0)
        # (3-element array specifying orientation of angular momentum vector)
        self.ang_mom = self.ang_mom_3D/ np.sum(_m)
        
        #
        # inclination orientation 
        #
        
        j=self.ang_mom/np.linalg.norm(self.ang_mom)
        #normalised specific angular momentum 
        
        x = np.array([1,2,3])
        x = x-(x.dot(j)*j) #make x orthogonal to j
        
        x/= np.linalg.norm(x) # normalise
        
        y = np.cross(j,x)#create 3rd vector - orth to x,j
        
        
        A = (x,y,j) # transformation matrix
        
        self.pgas_coo=np.dot(A,self.pgas_coo.T).T # change co-ordinates
        self.pgas_vel = np.dot(A,self.pgas_vel.T).T

        #
        # Apply same process to stellar particle type
        #
        
        self.pstar_coo=np.dot(A,self.pstar_coo.T).T  #change coordinates
        self.pstar_vel=np.dot(A,self.pstar_vel.T).T



    def rad_gen(self):
        self.gas_radial = np.sqrt((self.pgas_coo[:,0]**2)+(self.pgas_coo[:,1]**2))
        #print(self.gas_radial)
        self.star_radial = np.sqrt((self.pstar_coo[:,0]**2)+(self.pstar_coo[:,1]**2))
        #print(self.str_radial)



    def df_gen(self):
        #gas data -> dataframe -> dfg
        self.dfg = pd.DataFrame({"x":self.pgas_coo[:,0], "y": self.pgas_coo[:,1],"z":self.pgas_coo[:,2],"rad":self.gas_radial,"m":self.pgas_dens,"met":(12+np.log10(self.pgas_met))})
        #star data -> dataframe -> dfs
        self.dfs = pd.DataFrame({"x":self.pstar_coo[:,0], "y": self.pstar_coo[:,1],"z":self.pstar_coo[:,2],"rad":self.star_radial,"m":self.pstar_m,"met":(12+np.log10(self.pstar_met))})


    def rad_norm(self,factor):
        #normalise dataframe radial values -> all subhalos will have identical x scale set by factor input 
        self.dfg = self.dfg.round(3)
        self.dfg.rad = factor*((self.dfg.rad-self.dfg.rad.min())/(self.dfg.rad.max()-self.dfg.rad.min()))
        #
        self.dfs.rad = factor*((self.dfs.rad-self.dfs.rad.min())/(self.dfs.rad.max()-self.dfs.rad.min()))



    def fit_lin(self,dfin,pc):
        annuli = pc*self.crit_dist
        popt,pcov = curve_fit(linear_fit,dfin['rad'],dfin['met'])
        #apply curve_fit function to particle data, for either linear or curved fit, 
        dfin.sort_values(by='rad', inplace=True)
        dfin = dfin[dfin['rad']<annuli]
        medfit = medfilt(dfin['met'],kernel_size=21) 

        plt.figure(figsize=(15,10))
        #sort dataframe values by radial value -> plot clarity for line plotting
        plt.plot(dfin['rad'], medfit,'g-')
        #plt.plot(dfin['rad'], dfin['met'],'g+')
        plt.plot(dfin['rad'], linear_fit(dfin['rad'],*popt),'r--')
        plt.xlim(0,(10*pc))
        plt.ylim(8,12)
        plt.xlabel("Radial Distance (Normalised Code Units)")
        plt.ylabel("12+log10$(O/H)$")
        plt.title("Metallicity Gradient for {}({}-snap-{})".format(self.subID, self.simID, self.snapID))
        plt.tick_params(axis='both',which='both',direction='inout',length=15)
        filename = ("lin_fit_{}.png".format(self.subID))
        plt.savefig(filename)
        plt.close()
        
    def savgol(self,dfin,pc):
        annuli = pc*self.crit_dist
        #apply curve_fit function to particle data, for either linear or curved fit, 
        dfin.sort_values(by='rad', inplace=True)
        dfin = dfin[dfin['rad']<annuli]
        medfit = medfilt(dfin['met'],kernel_size=21) 
        
        interp_savgol = savgol_filter(dfin['met'], window_length=101, polyorder=3)

        plt.figure(figsize=(15,10))
        #sort dataframe values by radial value -> plot clarity for line plotting
        plt.plot(dfin['rad'], medfit,'b-')
        plt.plot(dfin['rad'], dfin['met'],'g+')
        plt.plot(dfin['rad'], interp_savgol,'r--')
        plt.xlim(0,(10*pc))
        plt.ylim(8,12)
        plt.xlabel("Radial Distance (Normalised Code Units)")
        plt.ylabel("12+log10$(O/H)$")
        plt.title("Metallicity Gradient for {}({}-snap-{})".format(self.subID, self.simID, self.snapID))
        plt.tick_params(axis='both',which='both',direction='inout',length=15)
        filename = ("savpng/lin_fit_{}.png".format(self.subID))
        plt.savefig(filename)
        plt.close()
    def bootleg_fit(self, runs):
        avals = []
        bvals = []
        pcov0 = []
        pcov1 = []
        for i in range(runs):
            df = self.dfg.sample(frac=0.1, replace=False)
            popt,pcov = curve_fit(linear_fit,df['rad'],df['met'])
            avals.append(popt[0])
            bvals.append(popt[1])
            pcov0.append(pcov[0])
            pcov1.append(pcov[1])
        return avals,bvals,pcov0, pcov1






#'''
df_in = pd.read_csv("rad.csv")

valid_id = df_in[df_in['sfr']>10e-1]
valid_id = valid_id[valid_id['radius']>9]
valid_id.to_csv("remove.csv")
valid_id = list(valid_id['ids'])
print(len(valid_id))
#'''

sub1 = galaxy("TNG50-1", 99, 8,'Y')
sub1.galcen()
sub1.ang_mom_align('gas')
sub1.rad_gen()
sub1.df_gen()
sub1.rad_norm(10)
sub1.savgol(sub1.dfg,0.75)
avals,bvals,pcov0,pcov1=sub1.bootleg_fit(20)
sub1.fit_lin(sub1.dfg, 1)
print("Min a {} : Max a: {} Mean a {}".format(min(avals),max(avals),np.mean(avals)))
print("Min b {} : Max b: {} Mean b {}".format(min(bvals),max(bvals),np.mean(bvals)))

'''

invalids =[]
def runlist(i):
    try:
        sub1 = galaxy("TNG50-1", 99, i,'Y')
        if sub1.test==0:
            print("sub invalid")
            invalids.append(i)
        else:
            sub1.galcen()
            sub1.ang_mom_align('gas')
            sub1.rad_gen()
            sub1.df_gen()
            sub1.rad_norm(10)
            sub1.savgol(sub1.dfg,0.75)
            print("Subhalo {} done".format(i))
    
    except OSError:
        print("corrupted file skipped")
        invalids.append(i)
    except KeyError:
        print("keyerror - skipped")
        invalids.append(i)

returns = Parallel(n_jobs=20)(delayed(runlist)(i) for i in valid_id)
filenames=[]
for i in invalids:
    filename = "fitpng/lin_fit_{}.png".format(i)
    filenames.append(filename)
for j in filenames:
    try:
        os.remove(j)
    except OSError:
        print('%%%woo')
'''

end = time.time()
print("Programme Runtime = {}".format(end-start))