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

class subhalo:
    def __init__(self,simID,snapID,subID,primeID):
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
        self.primesub = primeID
        #basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output' #path to simulation data on Noether
        basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'
        baseurl = 'http://www.tng-project.org/api/'+str(simID)+'/snapshots/'+str(snapID) # API url for simulation/snapshot information
        
        hubble =0.7
        
        #use api to easily read snapshot information
        redshift = 0 if self.snapID is 99 else 0.5 if self.snapID is 67 else 2 if snapID is 33 else print("invalid snapID: Not on Noether")
        
        #
        #redshift = 0 for snap99
        #redshift = 2.00202813925285 for snap 33
        #redshift = 0.503047523244883 for snap 67
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
        self.tot_met = all_fields['SubhaloGasMetallicity']
        self.m_tot = all_fields['SubhaloMass']
        self.totsfr = all_fields['SubhaloSFR']
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
        # Velocities  (N,3) km sqrt(scalefac)        
        # We convert these to pkpc (proper kpc), Msun and km/s, respectively
        crit_dist = 5 * self.Rhalf #30. # proper kpc
        self.crit_dist = crit_dist
        hcoldgas  = np.where( (gas['StarFormationRate'] > 0.0) & (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        self.test = len(hcoldgas)
        #print(hcoldgas)
        #print(len(hcoldgas1), len(hcoldgas2))
        #hcoldgas  = (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2)
        self.pgas_coo   = gas['Coordinates'][hcoldgas]/hubble / (1. + redshift)
        self.pgas_m     = gas['Masses'][hcoldgas] * 10**10 / hubble
        self.pgas_vel   = (gas['Velocities'][hcoldgas] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.conv_kms2kpcyr = (3.1558 / 3.08568) * 10**(-9)
        self.pgas_vel   = self.pgas_vel * self.conv_kms2kpcyr    #Convert to kpc/yr
        self.pgas_sfr   = gas['StarFormationRate'][hcoldgas]
        self.pgas_met   =gas['GFM_Metallicity'][hcoldgas]
        self.pgas_dens = gas['Density'][hcoldgas]
        self.pgas_sfr= gas['StarFormationRate'][hcoldgas]
        #print(all_fields.keys())
        # Load all stellar particle data
        stars = il.snapshot.loadSubhalo(basePath, snapID, subID, 'stars', fields=['Coordinates', 'Masses', 'Velocities','GFM_Metallicity' ])
        hstar = np.where( (np.sum((stars['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        self.pstar_coo   = stars['Coordinates'][hstar]/hubble / (1. + redshift)
        self.pstar_m     = stars['Masses'][hstar] * 10**10 / hubble
        self.pstar_vel   = (stars['Velocities'][hstar] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.pstar_vel   = self.pstar_vel * self.conv_kms2kpcyr
        self.pstar_met = stars['GFM_Metallicity'][hstar]
        
    def galcen_datagen(self):
        r'''
        Positions co-ordinates on frame of galactic centre of mass
        '''
        self.pgas_coo -= self.centre[None,:]
        self.pstar_coo -= self.centre[None,:]
        self.gas_radial = np.sqrt((self.pgas_coo[:,0]**2)+(self.pgas_coo[:,1]**2))
        
        df = pd.DataFrame({"x":self.pgas_coo[:,0],
                                   "y":self.pgas_coo[:,1],
                                   "z":self.pgas_coo[:,2],
                                   "rad": self.gas_radial,
                                   "mass":self.pgas_m,
                                   "dens":self.pgas_dens,
                                   "met":(self.pgas_met),
                                   "met2":(self.pgas_met ),
                                   "sfr":self.pgas_sfr})
        self.df = df
        return df
    
    def combfilter(self,dfin,scale):
        r'''
        Combined filter and normalisation to aid computational efficiency 

        '''
        df = dfin
        scaleheight = 0.1* self.stellarphotometricsrad
        df.rad = scale*((df.rad-df.rad.min())/(df.rad.max()-df.rad.min()))
        df = df[df['z']<scaleheight]
        self.df_norm = df
        return df
    
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
      

class dodirectory:
    def __init__(self,primeID):
        self.primeID = primeID
    def getlist(self):
        listpath = "files/historycutouts/evdir_{}/locals_{}.csv".format(self.primeID,self.primeID)
        df = pd.read_csv(listpath)
        snapshots = list(df['snapshots'])
        subhalos = list(df['subhalos']) 
        return snapshots,subhalos
    

def dosingle(sub,snap,prime):
    obj = subhalo('TNG50-1',snap,sub, prime)
    df = obj.galcen_datagen()
    df2 = obj.combfilter(df,2)
    subID = sub
    snapID = snap
    slope1,slope2 = obj.piecewise(df2,3)
    print("done for subhalo {} snapshot {}".format(sub,snap))
    return (subID, snapID,slope1,slope2)

def dodir(i):
    try:
        data = dodirectory(i)
        snapshots,subhalos = data.getlist()
        subs = [];snaps=[];s1=[];s2=[]
        for j in range(3):
            subID,snapID,slope1,slope2 = dosingle(subhalos[j],snapshots[j],i)
            subs.append(subID);snaps.append(snapID)
            s1.append(slope1);s2.append(slope2)
        df = pd.DataFrame({
            'subhalo':subs,
            'snapshot':snaps,
            'slope1':s1,
            'slope2':s2
        })
        #returns = Parallel(n_jobs=4)(delayed(dosingle)(subhalos[j],snapshots[j],i)for j in range(4))
        #df = pd.DataFrame(returns,columns = ['subhalo','snapshot','slope1','slope2'])
        fpath = "files/historycutouts/evdir_{}/localslope{}.csv".format(i,i)
        df.to_csv(fpath)
        return print("done for descendant {}".format(i))
    except OSError as e:
        return print(e)
    except TypeError as e:
        return print(e)
    except IndexError as e:
        return print(e)
    except ValueError as e:
        return print(e)
    

dfyay = pd.read_csv("traceids.csv")
ids = list(dfyay['id'])
returns = Parallel(n_jobs=20)(delayed(dodir)(i) for i in ids)