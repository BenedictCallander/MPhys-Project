# Visualisation.py 
# \-> contains condensed subhalo class to provide required analysis routines, with a number of visualisation plot routines provided
# Created:17/11/2022 
# Author: Benedict Callander 
# Respository https://github.com/btcallander/MPhys-Project (private)
#

#Plotting, numerical functions and dataset manipulation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pwlf

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


#set basic constants during initialisation for easy 
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
pd.options.mode.chained_assignment = None  # default='warn'




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
        #print(hcoldgas)
        #print(len(hcoldgas1), len(hcoldgas2))
        #hcoldgas  = (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2)
        self.test = len(hcoldgas)
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

    def galcen(self):
        r'''
        Positions co-ordinates on frame of galactic centre of mass
        '''
        self.pgas_coo -= self.centre[None,:]
        self.pstar_coo -= self.centre[None,:]
        
    def ang_mom_align(self, type):
        r'''
        Align subhalo with Z axis by calcualting angular velocity,momentum of particle cells and generating transformation matrix 
        
        INPUT: 
        
        Type: str
        
        Select either 'gas' or 'stars' for cell-type by which to perform alignment 
        '''
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
        r'''
        Transform Cartesian Co-ordinates of gas and star cells into radial cells -> create new dataset alongside existing XYZ
        '''
        self.gas_radial = np.sqrt((self.pgas_coo[:,0]**2)+(self.pgas_coo[:,1]**2))
        self.star_radial = np.sqrt((self.pstar_coo[:,0]**2)+(self.pstar_coo[:,1]))

    def df_gen(self,type,quant):
        r'''
        generate dataframe containing particle level data for subhalo object
        
        INPUTS:
        
        Type: str 
        
        String input to determine whether dataframe contains gas cell or star cell information 
        
        Quant: Str
        
        What parameters other than co-ordinates to include ('mass','dens','met',('comb' contains all ))
        
        '''
        
        #series of logical statements read two input parameters to generate dataframe suited for request type 
        if (type == 'gas'):
            if (quant == 'mass'):
                df = pd.DataFrame({"x":self.pgas_coo[:,0],
                                   "y":self.pgas_coo[:,1],
                                   "z":self.pgas_coo[:,2],
                                   "rad": self.gas_radial,
                                   "mass":self.pgas_m})
            elif (quant =='dens'):
                df = pd.DataFrame({"x":self.pgas_coo[:,0],
                                   "y":self.pgas_coo[:,1],
                                   "z":self.pgas_coo[:,2],
                                   "rad": self.gas_radial,
                                   "dens":self.pgas_dens})
            elif (quant =='met'):
                df = pd.DataFrame({"x":self.pgas_coo[:,0],
                                   "y":self.pgas_coo[:,1],
                                   "z":self.pgas_coo[:,2],
                                   "rad": self.gas_radial,
                                   "met":(self.pgas_met)})
            elif (quant =='comb'):
                df = pd.DataFrame({"x":self.pgas_coo[:,0],
                                   "y":self.pgas_coo[:,1],
                                   "z":self.pgas_coo[:,2],
                                   "rad": self.gas_radial,
                                   "mass":self.pgas_m,
                                   "dens":self.pgas_dens,
                                   "met":(self.pgas_met),
                                   "met2":(self.pgas_met ),
                                   "sfr":self.pgas_sfr})
        elif (type =='star'):
            if (quant == 'mass'):
                df = pd.DataFrame({"x":self.pstar_coo[:,0],
                                   "y":self.pstar_coo[:,1],
                                   "z": self.pstar_coo[:,2],
                                   "rad": self.star_radial,
                                   "mass": self.pstar_m})
            elif (quant =='met'):
                df = pd.DataFrame({"x":self.pstar_coo[:,0],
                                   "y":self.pstar_coo[:,1],
                                   "z": self.pstar_coo[:,2],
                                   "rad": self.star_radial,
                                   "met": self.pstar_met})
            elif (quant =='comb'):
                df = pd.DataFrame({"x":self.pstar_coo[:,0],
                                   "y":self.pstar_coo[:,1],
                                   "z": self.pstar_coo[:,2],
                                   "rad": self.star_radial,
                                   "mass": self.pstar_m,
                                   "met": self.pstar_met})
        
        self.df = df
        return df
    
    def gas_visual(self,dfin,lim):
        df = dfin
        df = df[abs(df['x'])<lim]
        df = df[abs(df['y'])<lim]
        df = df.groupby(['x','y'])['dens'].sum().reset_index()
        print(np.mean(df['dens']))
        print(df)
        plt.figure(figsize=(20,12))
        plt.style.use('dark_background')
        #plt.hist2d(df['x'],df['y'],weights =(df['dens']), bins=[2000,1500],cmap='magma',vmax=0.000006)#,vmin=0.0008,vmax=0.002)
        plt.scatter(df['x'],df['y'],c=10*(df['dens']),cmap='magma')#,vmax=0.0003)#,vmin=-7,vmax=-1)
        plt.xlabel('$\Delta x$ [kpc/h]')
        plt.ylabel('$\Delta y$ [kpc/h]')
        plt.colorbar(label='log10(Gas Mass)')
        plt.title('Gas Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
        #filename = 'Visuals/visuals/Mgass_{}_sub_{}.png'.format(self.simID, self.subID)
        filename = 'scatvis_{}.png'.format(self.subID)
        plt.savefig(filename)
        plt.close()
        
    def doublepiecewise(self,dfin,breakpoint1,breakpoint2):
        df = dfin#.sample(frac=0.01,replace=False)
        df.sort_values(by="rad",inplace = True)
        #med_data1 = medfilt((12+np.log10(df['met'])), kernel_size=11)
        x0 = np.array([min(df['rad']), breakpoint1,breakpoint2, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met']),weights=1/df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        slope1 = my_pwlf.slopes[0]
        slope2 = my_pwlf.slopes[1]
        slope3 = my_pwlf.slopes[2]
    
        #'''
        print("slopes are inner: {} middle:{} and outer:{}".format(slope1,slope2,slope3))
        xHat = np.linspace(min(df['rad']), max(df['rad']), num=10000)
        yHat = my_pwlf.predict(xHat)
        plt.figure(figsize=(20,12))
        plt.hist2d(df['rad'],12+np.log10(df['met']),bins=[200,200], weights=df['sfr'],cmap='binary')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        filename = '{}_sub_{}_doublebreak.png'.format(self.snapID, self.subID)
        plt.savefig(filename)
        plt.close()
        #'''
        return (slope1,slope2,slope3)

class csvvis:
    def __init__(self,fpath):
        self.path = fpath
    def slopegrad(self):
        df = self.df
        
        
'''
test subject= subhalo 275545
checkpoint 87768 on terminal
'''    

sub = subhalo('TNG50-1',99,96762)
sub.galcen()
sub.ang_mom_align('gas')
sub.rad_transform()
df =sub.df_gen('gas','comb')
print("current runtime {}".format(time.time()-start))
#sub.gas_visual(df,30)
sub.doublepiecewise(df,1,2)
#96762
#117250
end = time.time()
print("runtime = {}".format(end-start))

#np.log10