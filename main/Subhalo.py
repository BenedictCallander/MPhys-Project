#
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

#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il

#Own module containing utility functions 
import BCUTILS

# runtime calculation 
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
                                   "met2":(self.pgas_met/self.pgas_dens)})
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
        df = dfin
        scaleheight = 0.1* self.stellarphotometricsrad
        df.rad = scale*((df.rad-df.rad.min())/(df.rad.max()-df.rad.min()))
        df = df[df['z']<scaleheight]
        self.df_norm = df
        return df

    def AIC_test(self,dfin,breakpoint):
        r'''
        Function to determine whether linear or broken fit provides better representation of subhalo metallicity gradient 

        Calculate linear fit (popt,pcov and apply linear fit and radial data)

        save fit data -> calculate RSS between met value and 

        calculate split fit for first part (popt1, pcov1 (for first half))

        calculate split fit for second part (popt2,pcov2 (for second half))

        calculate RSS values and AIC values for each method, 
        conditional statement to determine which fit is better according to AIC 
        return statement/value which inidcates which fit is better -> also catalogue other data? 
        '''
        popt,pcov = curve_fit(UTILITY.linear_fit, dfin['rad'],dfin['met'])

        y1 = UTILITY.linear_fit(dfin['rad'],*popt)
        AIC = 1
        x0 = np.array([min(dfin['rad']), breakpoint, max(dfin['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(dfin['rad'], dfin['met'])
        my_pwlf.fit_with_breaks(x0)
        
        xHat = np.linspace(min(dfin['rad']), max(dfin['rad']), num=len(dfin['met']))
        yHat = my_pwlf.predict(xHat)
        RSS_linear = np.sum((dfin['met']-(y1))**2)
        RSS_broken = np.sum((dfin['met']-(yHat))**2)
        #print(RSS_broken)
        #print(RSS_linear)
        linear = 4+(len(y1)*np.log(RSS_linear))
        broken = 8+(len(yHat)*np.log(RSS_broken))
        linear = abs(linear)
        broken = abs(broken)
        if linear>broken:
            return 1 #broken provides better fit 
        elif linear<broken:
            return 2 #linear provides better fit 
        else:
            return 3

    def fit_linear(self,dfin, pc):
        r'''
        fit linear regression to subhalo metallicity gradient (y=m*x+c)

        INPUTS:

        Self:

        access to subhalo object properties 

        Dfin:

        Allows for specification of dataframe (if certain filter processes have been performed)

        PC:

        Max Radial distance of consideration (annuli)
        '''
        annuli = pc
        dfin = dfin[dfin['rad']<pc]  
        dfin.sort_values(by='rad',inplace = True)   
        popt,pcov = curve_fit(UTILITY.linear_fit, dfin['rad'],12+np.log10(dfin['met2']))
        med_data = medfilt(12+np.log10(dfin['met2']),kernel_size = 21)

        plt.figure(figsize=(20,12))
        plt.plot(dfin['rad'], med_data, 'r-')
        plt.plot(dfin['rad'], UTILITY.linear_fit(dfin['rad'],*popt))
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}O/H)$")
        filename = 'sub_{}_met_snap={}'.format(self.subID,self.snapID)
        plt.savefig(filename)
        plt.close()

    def fit_quad(self,dfin,pc):
        r'''
        Pseudocode
        Calculate quadratic fit (f(x)=a*x**2 + b*x + c)
        Take input of dataframe and fit quadratic trendline using scipy curve_fit 
        popt = [a,b,c]
        return -> plot png file saved to directory with ID, simulation and snapshot as filename/title

        '''
        annuli = pc
        dfin = dfin[dfin['rad']<pc]     
        popt,pcov = curve_fit(UTILITY.sq_fit, dfin['rad'],dfin['met'])
        med_data = medfilt(dfin['met'],kernel_size = 21)

        plt.figure(figsize=(20,12))
        plt.plot(dfin['rad'], med_data, 'r-')
        plt.plot(dfin['rad'], UTILITY.sq_fit(dfin['rad'],*popt))
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10} (\frac{O}{H})$")
        filename = 'quadraticfit/sub_{}_met_snap={}'.format(self.subID,self.snapID)
        plt.savefig(filename)
        plt.close()
    
    def broken_fit(self,dfin,breakpoint,pc):
        r'''
        Pseudocode
        Inputs:

        Dataframe (can be inherited from subhalo object)
        breakpoint (the point at which break in linear fit is placed)
        '''
        annuli = pc
        dfin.sort_values(by='rad',inplace = True)
        med_data1 = medfilt(12+np.log10(dfin['met2']), kernel_size=21)
        x0 = np.array([min(dfin['rad']), breakpoint, max(dfin['rad'])])
        break1 = list(dfin[dfin['rad']<breakpoint])
        break2 = list(dfin[dfin['rad']>=breakpoint])
        my_pwlf = pwlf.PiecewiseLinFit(dfin['rad'], 12+np.log10(dfin['met2']))
        my_pwlf.fit_with_breaks(x0)
        
        xHat = np.linspace(min(dfin['rad']), max(dfin['rad']), num=10000)
        yHat = my_pwlf.predict(xHat)
        
        plt.figure(figsize=(20,12))
        plt.plot(dfin['rad'], med_data1, 'b--')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        filename = 'sub_{}_break_snapshot_{}.png'.format(self.subID, self.snapID)
        plt.savefig(filename)
        plt.close()
    
    def savgol_smooth(self,dfin,pc,scatter_q):
        r'''
        INPUTS:
        
        dfin: pandas DataFrame (must contain 'rad', 'met' columns)

        scatter_q: str
        
        'Y' or 'N' : determines whether scatterplot also included (demonstration purposes)

        Calculates a Savitzky-Golay filter -> used to smooth metallicity data for clearer representation of the gradient 

        WARN: values calculated by the Savitzky-Golay filter CANNOT be used to calculate scientific results 

        returns:  plot png file saved to directory with ID, simID and snapID as filename/title

        Optional Extras: could include conditions to plot additional information onto graph such as linear/broken/quadratic fits? 
        '''
        annuli = pc
        dfin.sort_values(by='rad',inplace=True)
        dfin=dfin[dfin['rad']<annuli]
        interp_savgol = savgol_filter(dfin['met'],window_length=101,polyorder=3)
        plt.figure(figsize=(15,10))
        if scatter_q=='Y':
            plt.plot(dfin['rad'],dfin['met'],c='g',marker = '+', markersize=2)
        elif scatter_q=='N':
            print("No scatter plot")
        plt.plot(dfin['rad'], interp_savgol, 'r--')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        filename = 'Visuals/savgolfits/sub_{}_savgol_snap_{}'.format(self.subID, self.snapID)
        plt.savefig(filename)
        plt.close()

    def bootstrap_test(self,dfin,runs,frac):
        r'''
        statistical test to test impact of outliers on linear fit

        INPUTS:
        
        dfin: pandas DataFrame
        
        dataframe containing metallicity and radius data for subhalo 

        runs: int
        
        number of bootstrap tests to be run 

        frac: 0<frac<=1

        Fraction of datapoints to be sampled in each run 

        RETURNS:
        m,b, pcov[0],[1]
        gradient, intercept and goodness of fit parameters from curve_fit 

        '''
        avals = []; bvals = [] ; pcov0 = []; pcov1=[]
        df = dfin
        for i in range(runs):
            sample = df.sample(frac=frac, replace=False)
            popt,pcov = curve_fit(UTILITY.linear_fit,sample['rad'],sample['met'])
            avals.append(popt[0])
            bvals.append(popt[1])
            pcov0.append(pcov[0])
            pcov1.append(pcov[1])
        return avals,bvals,pcov0, pcov1
 
    def gas_visualisation(self, dfin, decp):
        r'''
        Inputs -> dataframe containing gas density data (aligned to z axis)

        flattens gas density data by decp (0 = 1kpc box, 1 .1kpc , 2 = .01kpc etc)
        
        plot visual by scatter of density values -> see other codes for plot syntax
        
        '''
        df = dfin.round(decp)
        df = df.groupby(['x','y'])['dens'].sum().reset_index()
        plt.figure(figsize=(20,12), dpi=500)
        plt.style.use('dark_background')
        plt.hist2d(df['x'],df['y'],weights = np.log10(df['dens']), bins=[1000,1000],cmap='magma')
        #plt.scatter(df['x'],-df['y'],c=(np.log10(df['m'])),cmap='inferno', vmin=(min(np.log10(df['m']))),vmax =(max(np.log10(df['m']))))
        plt.xlabel('$\Delta x$ [kpc/h]')
        plt.ylabel('$\Delta y$ [kpc/h]')
        #plt.colorbar(label='log10(Gas Mass)')
        plt.title('Gas Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
        #filename = 'Visuals/visuals/Mgass_{}_sub_{}.png'.format(self.simID, self.subID)
        filename = 'histtest.png'
        plt.savefig(filename)
        plt.close()
    
    def slopegen(self,breakpoint,dfin):
        r'''
        slopegen: 
        
        Determine best fit (linear,broken) for subhalo Metallicity Gradient, and return gradient (inner gradient returned in case of broken fit)

        INPUTS:
        
        self: self
        
        access data contained within subhalo object
        
        breakpoint: position of breakpoint in piecewise linear fit (constant)
        
        save fit data -> calculate RSS between met value and 

        calculate split fit for first part (popt1, pcov1 (for first half))

        calculate split fit for second part (popt2,pcov2 (for second half))

        calculate RSS values and AIC values for each method, 
        conditional statement to determine which fit is better according to AIC 
        return statement/value which inidcates which fit is better -> also catalogue other data? 
        '''
        df=dfin
        popt,pcov = curve_fit(UTILITY.linear_fit, df['rad'],df['met'])

        y1 = UTILITY.linear_fit(df['rad'],*popt)
        AIC = 1
        x0 = np.array([min(df['rad']), breakpoint, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], df['met'],disp_res = True)
        my_pwlf.fit_with_breaks(x0)
        brokenslope = my_pwlf.slopes[0]
        
        xHat = np.linspace(min(df['rad']), max(df['rad']), num=len(df['met']))
        yHat = my_pwlf.predict(xHat)
        RSS_linear = np.sum((df['met']-(y1))**2)
        RSS_broken = np.sum((df['met']-(yHat))**2)
        #print(RSS_broken)
        #print(RSS_linear)
        linear = 4+(len(y1)*np.log(RSS_linear))
        broken = 8+(len(yHat)*np.log(RSS_broken))
        linear = abs(linear)
        broken = abs(broken)
        
        met = self.tot_met
        mass = self.m_tot
        sfr = self.totsfr
        
        if linear>broken:
            return (brokenslope,met,mass,sfr)
        elif broken>linear:
            return (popt[0],met,mass,sfr)
        else:
            return (popt[0],met,mass,sfr)

#-------------------------------------------------------------------------------------------------------------------------------------|
#set simulation snapshot and get IDS for star forming subhalos -> pass to function to create object for each/ perform analysis on all |
#-------------------------------------------------------------------------------------------------------------------------------------|

sim = 99
dfin = BCUTILS.subhalo_classification(sim)
valid_id = list(dfin['id'])

#--------------------------------------------------------------------------------------------------------------------------------------|
#Function to generate subhalo object using class and specify analysis functions to be run -> function allows parallelisation to be used|
#--------------------------------------------------------------------------------------------------------------------------------------|

def subhalo_analysis(i):
    try:
        sub = subhalo("TNG50-1",99,i)
        if sub.test<4:
            print("not enough gas cells to continue")
        else:
            sub.galcen()
            sub.ang_mom_align('gas')
            sub.rad_transform()
            dfg = sub.df_gen('gas','comb')
            dfg2 = sub.combfilter(dfg,10)
            idval = i
            slope,met,mass,sfr = sub.slopegen(5,dfg2)
            print("subhalo {} calculated: current runtime time: {}".format(i,(time.time()-start)))
            return (slope,met,mass,idval,sfr)
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
        return print(e)
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

returns = Parallel(n_jobs= 60)(delayed(subhalo_analysis)(i) for i in valid_id)
df2=pd.DataFrame(returns,columns=['slope','met','mass(wrongunits)','id','sfr'])
df2.insert(5,'massplot', dfin['mass'],True)
df2.to_csv("csv/tng99slopes.csv")

#------------------------------------------------------------------------------------------------------------------------------|
# Pass dataframes into BCUTILS MSfilter function to create dataset containing only main sequence subhalos for separate analysis|
#------------------------------------------------------------------------------------------------------------------------------|
BCUTILS.MSfilter(dfin,df2,'csv/tng99MAIN.csv')

end = time.time()
print('runtime = {}s'.format(end-start))