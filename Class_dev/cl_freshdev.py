#
# This wasnt a good idea 
#

import requests
import h5py 
import numpy as np 
import matplotlib.pyplot as plt 
import scipy
import pandas as pd 
import time 
from joblib import Parallel, delayed
from scipy.signal import medfilt
from scipy.optimize import curve_fit 

headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
baseurl='https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/'

def get(path, params = None):
    '''
    API request function - expansion upon requests.get() to provide error codes and .hdf5 file reading 
    Inputs:
        - Path - URL for request to fetch data from 
        - Params - additional parameters that can be used (i.e. search and limit in the case of the ILLUSTRIS TNG API  )
    '''
    #Make API request for path (needs to be str variable) including header (API KEY) and other params
    r = requests.get(path, params=params,headers=headers) 

    #check for HTTP error codes
    r.raise_for_status()

    #detect content type for function to read 

    if r.headers['content-type'] == 'application/json':
        return r.json()
    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r



class APIutils:
    def __init__(self, url,):
        self.url=url
    def masssfrplot(self,criteria, filename):
        baseurl = self.url
        url = baseurl+str(criteria)
        subsdata=get(url)
        ids = []
        mass=[]
        sfr=[]
        valid_ids = [ subsdata['results'][i]['id'] for i in range(17553)]
        for i,id in enumerate(valid_ids):
            mass.append(subsdata['results'][i]['mass_log_msun'])
            ids.append(subsdata['results'][i]['id'])
            sfr.append(subsdata['results'][i]['sfr'])

        plt.figure(figsize=(15,10))
        plt.plot(mass, np.log10(sfr), 'k+')
        plt.xlabel("Mass (Log Msun)")
        plt.ylabel("log(SFR)")
        plt.savefig(filename)
        plt.close()


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