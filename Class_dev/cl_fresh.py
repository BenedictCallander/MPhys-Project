#
# Fresh Development -> eliminate bugs and restructure programme 
#
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

headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
#basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output'
pd.options.mode.chained_assignment = None  # default='warn'



class UTILITY:
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
    def __init__(self,simID,snapID,subID):
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
        #hcoldgas  = np.where( (gas['StarFormationRate'] > 0.) & (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
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

    def rad_transform(self):
        self.gas_radial = np.sqrt((self.pgas_coo[:,0]**2)+(self.pgas_coo[:,1]**2))
        self.star_radial = np.sqrt((self.pstar_coo[:,0]**2)+(self.pstar_coo[:,1]))

    def df_gen(self,type,quant):
        if (type == 'gas'):
            if (quant == 'mass'):
                df = pd.DataFrame({"x":self.pgas_coo[:,0],"y":self.pgas_coo[:,1], "z":self.pgas_coo[:,2], "rad": self.gas_radial, "mass":self.pgas_m})
            elif (quant =='dens'):
                df = pd.DataFrame({"x":self.pgas_coo[:,0],"y":self.pgas_coo[:,1], "z":self.pgas_coo[:,2], "rad": self.gas_radial, "dens":self.pgas_dens})
            elif (quant =='met'):
                df = pd.DataFrame({"x":self.pgas_coo[:,0],"y":self.pgas_coo[:,1], "z":self.pgas_coo[:,2], "rad": self.gas_radial, "met":self.pgas_met})
            elif (quant =='comb'):
                df = pd.DataFrame({"x":self.pgas_coo[:,0],"y":self.pgas_coo[:,1], "z":self.pgas_coo[:,2], "rad": self.gas_radial, "mass":self.pgas_m, "dens":self.pgas_dens, "met": self.pgas_met})
        elif (type =='star'):
            if (quant == 'mass'):
                df = pd.DataFrame({"x":self.pstar_coo[:,0], "y":self.pstar_coo[:,1], "z": self.pstar_coo[:,2], "rad": self.star_radial, "mass": self.pstar_m})
            elif (quant =='met'):
                df = pd.DataFrame({"x":self.pstar_coo[:,0], "y":self.pstar_coo[:,1], "z": self.pstar_coo[:,2], "rad": self.star_radial, "met": self.pstar_met})
            elif (quant =='comb'):
                df = pd.DataFrame({"x":self.pstar_coo[:,0], "y":self.pstar_coo[:,1], "z": self.pstar_coo[:,2], "rad": self.star_radial, "mass": self.pstar_m, "met": self.pstar_met})
        self.df = df
        return df

    def rad_norm(self, dfin, scale):
        # INPUTS
        #dfin -> dataframe with keyword 'rad' to be normalised to code units
        #scale -> scale of code units -> size of axis (normalisation completes to between 0 and 1)
        df = dfin
        df.rad = scale*((df.rad-df.rad.min())/(df.rad.max()-df.rad.min()))
        self.df_norm = df
        return df
    
    def z_filter(self, dfin):
        #Takes Stellar photometrics radius and relation of 0.1* to filter height along z axis 
        scaleheight = 0.1* self.stellarphotometricsrad
        dfout = dfin[dfin['z']<scaleheight]

        return dfout

    def AIC_test(self):
        '''
        Pseudocode 

        Calculate linear fit (popt,pcov and apply linear fit and radial data)

        save fit data -> calculate RSS between met value and 

        calculate split fit for first part (popt1, pcov1 (for first half))

        calculate split fit for second part (popt2,pcov2 (for second half))

        calculate RSS values and AIC values for each method, 
        conditional statement to determine which fit is better according to AIC 
        return statement/value which inidcates which fit is better -> also catalogue other data? 
        '''
        
    def fit_linear(self,dfin):
        '''
        Pseudocode
        Calculate linear fit (f(x) = a*x+b)
        Take input of dataframe and fit linear trendline using scipy curve_fit 
        (popt,pcov) (popt[0] = gradient) , (popt[1]=intercept)

        return -> plot png file saved to directory with ID, simulation and snapshot as filename/title
        '''
    
    def fit_quad(self,dfin):
        '''
        Pseudocode
        Calculate quadratic fit (f(x)=a*x**2 + b*x + c)
        Take input of dataframe and fit quadratic trendline using scipy curve_fit 
        popt = [a,b,c]


        return -> plot png file saved to directory with ID, simulation and snapshot as filename/title

        '''
    
    def broken_fit(self,dfin,breakpoint):
        '''
        Pseudocode
        Inputs:
        Dataframe (generated in galaxy class )
        breakpoint (the point at which break in linear fit is placed)

        process
        split DF into 2 
        radfilt 1 = dfin[dfin['rad']<breakpoint]
        radfilt 2 = dfin[dfin['rad']>breakpoint]

        calculate popt, pcov linear fit for both datapoints 

        save popt, pcov values ?

        plot data with broken fit overlaid (use median filter for metallicity data?)
        
        '''
    
    def savgol_smooth(self,dfin):
        '''
        Pseudocode

        Takes input of dataframe with 'rad' and 'met' keywords, which represent radius and metlalicity 
        Calculates a Savitzky-Golay filter -> used to smooth metallicity data for clearer representation of the gradient 

        WARN - values calculated by the Savitzky-Golay filter cannot be used to calculate scientific results 

        return -> plot png file saved to directory with ID, simID and snapID as filename/title

        Optional Extras -> could include conditions to plot additional information onto graph such as linear/broken/quadratic fits? 
        '''
    
    def bootstrap_test(self,dfin,type):
        '''
        Pseudocode 

        Takes input of dataframe and applies bootstrapping depending on type (either linear or quadratic)

        samples dataframe keeping only fraction of values (random) each time -> 

        compare means, range 

        '''
    
    def gas_visualisation(self, dfin, decp):
        '''
        Psuedocode
        Inputs -> dataframe containing gas density data (aligned to z axis)

        flattens gas density data by decp (0 = 1kpc box, 1 .1kpc , 2 = .01kpc etc)
        
        plot visual by scatter of density values -> see other codes for plot syntax
        
        '''
    
