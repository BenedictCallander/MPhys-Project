#Plotting, numerical functions and dataset manipulation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pwlf
import scipy.stats as stats
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
start = time.time()

basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'

fields = ['SubhaloMassType','SubhaloSFR']
subhalos = il.groupcat.loadSubhalos(basePath,99,fields=fields)
mass_msun = subhalos['SubhaloMassType'][:,4] * 1e10/0.7
sfr = subhalos['SubhaloSFR']


class subhalo:
    def __init__(self,simID,snapID,subID):
        r'''
        Initialise subhalo object
        
        Inputs:
        
        simID: str
        
        Illustris TNG Simulation ID (e.g. 'TNG50-1', 'TNG100-1')
        
        snapID: float, int
        
        ID number of snapshot containing subhalo of study 
        
        subID: float, int
        
        ID number of target subhalo  
        '''
        np.seterr(divide='ignore', invalid='ignore')
        #object creation requires only 3 input parameters for galaxy selection 
        self.simID = simID #simulation used (TNG100/TNG50)
        self.snapID = snapID #snapshot of study - can be either number of z= (z= required in parenthesis )
        self.subID = subID #subhalo ID - could be iterated 
        #basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output' #path to simulation data on Noether
        basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'
        baseurl = 'http://www.tng-project.org/api/'+str(simID)+'/snapshots/'+str(snapID) # API url for simulation/snapshot information
        
        hubble =0.7
        
        #use api to easily read snapshot information
        redshift = 0 if self.snapID is 99 else 0.5 if self.snapID is 67 else 2 if snapID is 33 else print("invalid snapID: Not on Noether")
        
        self.redshift =redshift  # store redshift value as attribute
        
        scalefac = 1./(1.+redshift) #calculate scale factor
        
        #
        # Read Subhalo level info 
        #
        
        ptNumGas = il.snapshot.partTypeNum('gas') #determine index designation for each particle type
        ptNumStars = il.snapshot.partTypeNum('stars')
        #
        #pull all data for specific subhalo 
        #
        all_fields= il.groupcat.loadSingle(basePath, snapID, subhaloID = subID)
        self.test=all_fields['SubhaloMassInRadType'][ptNumGas]
        self.tot_met = all_fields['SubhaloGasMetallicity']
        self.tot_met = 8.69 + np.log10(self.tot_met) - np.log10(0.0127)
        self.rad_met = all_fields['SubhaloGasMetallicityHalfRad']
        self.rad_met = 8.69 + np.log10(self.rad_met) - np.log10(0.0127)
        self.m_tot = all_fields['SubhaloMass']
        self.totsfr = all_fields['SubhaloSFR']
        self.radsfr = all_fields['SubhaloSFRinHalfRad']
        self.totMstar = all_fields['SubhaloMassType'][ptNumStars]
        self.Mgas  = all_fields['SubhaloMassInRadType'][ptNumGas]
        self.Mstar = all_fields['SubhaloMassInRadType'][ptNumStars]
        # Coordinate of particle with minimum binding energy (converted from ckpc/h to kpc)
        self.centre = all_fields['SubhaloPos']/hubble / (1. + redshift)  # 3-element array [units: proper kpc]
        # Adopt the 3D half-stellar-mass radius
        self.Rhalf  = all_fields['SubhaloHalfmassRadType'][ptNumStars]/hubble / (1. + redshift)  
        self.stellarphotometricsrad = all_fields['SubhaloStellarPhotometricsRad']
        self.Rhalfdouble  = all_fields['SubhaloHalfmassRadType']
        # [units: proper kpc] (quantified in 3D)
        # Load all the relevant particle level info
        
def dosingle(i):
    sub = subhalo('TNG50-1',33,i)
    totmet = sub.tot_met
    radmet=sub.rad_met
    Mstar = sub.totMstar
    MstarRad = sub.Mstar
    totsfr = sub.totsfr
    radsfr = sub.radsfr
    subID = i
    print("done for {}".format(i))
    return (totmet,radmet,Mstar,MstarRad,totsfr,radsfr,subID)


dfin = pd.read_csv("tng33.csv")
ids = list(dfin['id'])
returns = Parallel(n_jobs=50)(delayed(dosingle)(i)for i in ids)
df2 = pd.DataFrame(returns, columns = ['met','metrad','mass','massrad','sfr','sfrrad','id'])
df2.to_csv("tng33subhalos.csv")
print("runtime = {}".format(time.time()-start))