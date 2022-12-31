# Subhalo.py 
# \-> contains subhalo class and subsequent analysis function (which runs all desired operations on each subhalo object)
# Created:17/11/2022 
# Author: Benedict Callander 
# Respository https://github.com/btcallander/MPhys-Project (private)
#

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


#set basic constants during initialisation for easy 
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
pd.options.mode.chained_assignment = None  # default='warn'

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
            filename = r.headers['content-disposition'].split("filename=")[1]
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

#
# GALAXY CLASS -> contains all subhalo analysis functions 
#

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
        # [units: proper kpc] (quantified in 3D)
        # Load all the relevant particle level info
        gas = il.snapshot.loadSubhalo(basePath, snapID, subID, 'gas', fields=['Coordinates', 'Masses','Density','Velocities', 'StarFormationRate','GFM_Metallicity','GFM_Metals'])
        # dimensions and units (see https://www.tng-project.org/data/docs/specifications/#parttype0):
        # Coordinates (N,3) ckpc/h   where ckps stands for co-moving kpc
        # Masses      (N)   10**10 Msun/h
        # Velocities  (N,3) km sqrt(scalefac)        
        # We convert these to pkpc (proper kpc), Msun and km/s, respectively
        crit_dist = 5 * self.Rhalf #30. # proper kpc
        self.crit_dist = crit_dist
        hcoldgas  = np.where((gas['StarFormationRate'] > 0.0) & (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        #hcoldgas  = np.where((np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        self.test = len(hcoldgas)
        #print(hcoldgas)
        #print(len(hcoldgas1), len(hcoldgas2))
        self.pgas_coo   = gas['Coordinates'][hcoldgas]/hubble / (1. + redshift)
        self.pgas_m     = gas['Masses'][hcoldgas] * 10**10 / hubble
        self.pgas_vel   = (gas['Velocities'][hcoldgas] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.conv_kms2kpcyr = (3.1558 / 3.08568) * 10**(-9)
        self.pgas_vel   = self.pgas_vel * self.conv_kms2kpcyr    #Convert to kpc/yr
        self.pgas_sfr   = gas['StarFormationRate'][hcoldgas]
        self.pgas_met   =gas['GFM_Metallicity'][hcoldgas]
        self.pgas_met = 8.69 + np.log10(self.pgas_met) - np.log10(0.0127)

        self.pgas_dens = gas['Density'][hcoldgas]
        self.pgas_sfr= gas['StarFormationRate'][hcoldgas]
        #self.abundance_O =gas['GFM_Metals'][:,4][hcoldgas]
        #self.abundance_H = gas['GFM_Metals'][:,1][hcoldgas]
        #self.pgas_met= self.abundance_O/self.abundance_H
        # 
        #Load all stellar particle data
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

    def rad_norm(self, dfin, scale):
        r'''
        Normalise radial co-ordinates to between (0->1)
        
        INPUTS:
        
        dfin: pandas DataFrame
        
        dataframe (generated by subhalo class) containing particle level data for subhalo object
        
        scale: int
        
        Scale of radial co-ordinates (transforms normalisation by linear factor)
        
        '''
        # INPUTS
        #dfin -> dataframe with keyword 'rad' to be normalised to code units
        #scale -> scale of code units -> size of axis (normalisation completes to between 0 and 1)
        df = dfin
        df.rad = scale*((df.rad-df.rad.min())/(df.rad.max()-df.rad.min()))
        self.df_norm = df
        return df
    
    def z_filter(self, dfin):
        r'''
        filter out cells at high z distances from galactic plane which do not contain significant gas to reduce number of cells computed
        '''
        #Takes Stellar photometrics radius and relation of 0.1* to filter height along z axis 
        scaleheight = 0.1* self.stellarphotometricsrad
        dfout = dfin[dfin['z']<scaleheight]

        return dfout
    
    def combfilter(self,dfin,scale):
        r'''
        Combined filter and normalisation to aid computational efficiency 

        '''
        df = dfin.copy()
        self.breakpoint = 10*((self.Rhalf-df.rad.min())/(df.rad.max()-df.rad.min()))
        scaleheight = 0.1* self.stellarphotometricsrad
        df = df[df['z']<scaleheight]
        df.rad = scale*((df.rad-df.rad.min())/(df.rad.max()-df.rad.min()))
        self.df_norm = df
        return df

    def median_profile(self,dfin):
        df = dfin.copy()
        df.sort_values(by='rad',inplace=True)
        ranges = np.arange(df.rad.min() - 0.1, df.rad.max() + 0.1, 0.1)  
        groups = df.groupby(pd.cut(df.rad, ranges))
        medianout = groups.median()
        medianout = medianout.dropna()
        return medianout

    def broken_fit(self,dfin):
        df = dfin.copy()
        
        x0 = np.array([df.rad.min(),self.breakpoint,df.rad.max()])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], (df['met']),weights=1/df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        slope1 = my_pwlf.slopes[0]
        slope2 = my_pwlf.slopes[1]
        return(slope1,slope2)
    
    def linear_fit(self,dfin):
        df = dfin.copy()
        popt,pcov = curve_fit(UTILITY.linear_fit, df['rad'],(df['met']),sigma = 1/df['sfr'], absolute_sigma= True)
           
        return (popt[0],df.rad.max())
    
    def getbreaks(self,dfin):
        df = dfin.copy()
        x = list(df['rad'])
        y = list(df['met'])
        piecewise_model = pwlf.PiecewiseLinFit(x, y)
        num_breaks = 1
        piecewise_model.fit(num_breaks)
        breaks = piecewise_model.fit_guess([self.breakpoint])
        self.breakpoint2 = breaks[1]
        print(self.breakpoint2)
        return breaks[1]
    
    def fit_and_compare(self,dfin):
    # fit linear model
        df = dfin.copy()
        x = list(df['rad'])
        y = list(df['met'])
        linear_model = np.polyfit(x, y, 1)
        linear_fit = np.poly1d(linear_model)

        # fit piecewise linear model
        piecewise_model = pwlf.PiecewiseLinFit(x, y)
        num_breaks = 1
        piecewise_model.fit(num_breaks)
        breaks = piecewise_model.fit_guess([6.0])
        piecewise_fit = piecewise_model.predict

        # compare models using AIC
        linear_residuals = y - linear_fit(x)
        linear_sse = np.sum(linear_residuals**2)
        linear_k = 2 # number of parameters in linear model (slope and intercept)
        linear_n = len(x) # number of data points
        linear_aic = 2*linear_k + linear_n * np.log(linear_sse/linear_n)

        piecewise_residuals = y - piecewise_fit(x)
        piecewise_sse = np.sum(piecewise_residuals**2)
        piecewise_k = num_breaks + 1 # number of parameters in piecewise model (slope and intercept for each line segment)
        piecewise_aic = 2*piecewise_k + linear_n * np.log(piecewise_sse/linear_n)

        if linear_aic < piecewise_aic:
            #print("Linear fit is better according to AIC")
            return 1
        else:
            #print("Piecewise fit is better according to AIC")
            return 2

    def bootstrap(self,dfin,runs,frac):
        df = dfin.copy()
        self.gradients = []
        for i in range(runs):
            sample = df.sample(frac=frac, replace=False)
            popt,pcov = curve_fit(UTILITY.linear_fit, sample['rad'],(sample['met']),sigma = 1/sample['sfr'], absolute_sigma= True)
            self.gradients.append(popt[0])
        
        minval = min(self.gradients)
        maxval = max(self.gradients)
        
        return minval,maxval            
#-------------------------------------------------------------------------------------------------------------------------------------|
#set simulation snapshot and get IDS for star forming subhalos -> pass to function to create object for each/ perform analysis on all |
#-------------------------------------------------------------------------------------------------------------------------------------|

#--------------------------------------------------------------------------------------------------------------------------------------|
#Function to generate subhalo object using class and specify analysis functions to be run -> function allows parallelisation to be used|
#--------------------------------------------------------------------------------------------------------------------------------------|
def subhalo_analysis(i):
    try:
        sub = subhalo("TNG50-1",33,i)
        if sub.test<10:
            print("not enough gas cells to continue")
        else:
            sub.galcen()
            sub.ang_mom_align('gas')
            sub.rad_transform()
            dfg = sub.df_gen('gas','comb')
            dfg2 = sub.combfilter(dfg,10)
            dfg3 = sub.median_profile(dfg2)
            slope, maxrad = sub.linear_fit(dfg2)
            slope1,slope2 = sub.broken_fit(dfg3)
            #minval,maxval=sub.bootstrap(dfg2,20,0.1)
            #change=100*((slope-minval)/(maxval-minval))
            #print("max-min = {}".format(change))
            fname = ("332slopes/slopedata_{}.csv".format(i))
            dfg3.to_csv(fname)
            idval=i
            met = sub.tot_met
            sfr = sub.totsfr
            Rhalf = sub.Rhalf
            mass = sub.Mstar
            print("subhalo {} calculated: current runtime: {}".format(i,(time.time()-start)))
            #return (met,idval,sfr,slope1,slope2,Rhalf)
            return (mass,met,idval,sfr,slope,Rhalf,slope1,slope2)

    except ValueError as e:
        #fname = "errors/errors{}.txt".format(i)
        #f = open(fname,"w")
        #f.write("errorcode: {} for subhalo {} \n".format(str(e),i))
        #f.close
        return print(e)
    
    except KeyError as e:
        #fname = "errors/errors{}.txt".format(i)
        #f = open(fname,"w")
        #f.write("errorcode: {} for subhalo {} \n".format(str(e),i))
        #f.close        
        return print(e,i)
    
    except OSError as e:
        #fname = "errors/errors{}.txt".format(i)
        #f = open(fname,"w")
        #f.write("errorcode: {} for subhalo {} \n".format(str(e),i))
        #f.close         
        return print(e)
    
    except TypeError as e:
        #fname = "errors/errors{}.txt".format(i)
        #f = open(fname,"w")
        #f.write("errorcode: {} for subhalo {} \n".format(str(e),i))
        #f.close         
        return print(e)

#-----------------------------------------------------------------------------------------------------------------------------------------|
#Call function in paralell computation to simultaneously perform analysis on (n_jobs) subhalos, write desired properties to dataframe->csv|
#-----------------------------------------------------------------------------------------------------------------------------------------|
'''
sub = subhalo("TNG50-1",99,37)
sub.galcen()
sub.ang_mom_align('gas')
sub.rad_transform()
dfg = sub.df_gen('gas','comb')
dfg2 = sub.combfilter(dfg,10)
dfg3 = sub.median_profile(dfg2)
sub.broken_fit(dfg3)
'''

#dfg3.to_csv("sfr2.csv")
sim = 99
dfin = pd.read_csv("tng33MS.csv")
dfin = dfin.dropna()
#pd.read_csv("csv/tng33MAIN.csv")
valid_id = list(dfin['id'])
returns = Parallel(n_jobs= 40)(delayed(subhalo_analysis)(i) for i in valid_id)
#df2=pd.DataFrame(returns,columns=['met','id','sfr','slope','maxrad','minval','maxval','change','Rhalf'])
df2=pd.DataFrame(returns,columns=['mass','met','id','sfr','slope','Rhalf','slope1','slope2'])
#df2 =pd.Dataframe(returns,columns = ['met','idval','sfr','minval','maxval','Rhalf'])
#df2.insert(5,'mass', dfin['mass'],True)
df2.to_csv("csv/tng33slopes.csv")




#------------------------------------------------------------------------------------------------------------------------------|
# Pass dataframes into BCUTILS MSfilter function to create dataset containing only main sequence subhalos for separate analysis|
#------------------------------------------------------------------------------------------------------------------------------|
end = time.time()
print('runtime = {}s'.format(end-start))

'''
Main sequence number of subhalos: 
snap33:8162
snap67:1581
snap99:1409

Main Sequence Line (y=10**(2x-20.5))

Number of star forming subhalos: (total subhalos processed)

snap33:77655
snap67:26980
snap99:17553

Author 	    G(dex/kpc) 	    A8.5         	c
Shaver 	    -0.05 ± 0.01 	8.77 ± 0.14 	-0.69
Simpson 	-0.06 ± 0.02 	8.57 ± 0.15 	-0.56
Afflerbach 	-0.06 ± 0.01 	8.61 ± 0.09 	-0.66
Maciel   	-0.07 ± 0.01 	8.66 ± 0.06 	-0.70
Fesen 	    -0.04 ± 0.03 	8.63 ± 0.32 	-0.31
Composite 	-0.06 ± 0.01

'''


'''
Keys for overview plots:

SubhaloGasMetallicity -  Mass-weighted average metallicity (Mz/Mtot, 
                        where Z = any element above He) of the gas cells bound to this Subhalo,
                        but restricted to cells within twice the stellar half mass radius. 
                        
SubhaloMass -  Total mass of all member particle/cells which are bound to this Subhalo, of all types.
                Particle/cells bound to subhaloes of this Subhalo are NOT accounted for. 


Assessors:

Mathlin: Cares about the physics - must demonstrate complete and comprehensive knowledge of background physics
        - what is the key physics in the intro 
        
        - clartity - GM not practicing astrophysicist so think of clarity and linguo
'''