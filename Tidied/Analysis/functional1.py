#
# Fresh Development -> eliminate bugs and restructure programme 
#

import time  # runtime calculation import numpy as np #data handling
from random import random
from re import sub  # http logging for debugging purpouses
import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests  # obtain data from API server
import pwlf

from joblib import Parallel, delayed
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter

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

#
# GALAXY CLASS -> contains all subhalo analysis functions 
#

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
                                   "met":12+np.log10(self.pgas_met)})
            elif (quant =='comb'):
                df = pd.DataFrame({"x":self.pgas_coo[:,0],
                                   "y":self.pgas_coo[:,1],
                                   "z":self.pgas_coo[:,2],
                                   "rad": self.gas_radial,
                                   "mass":self.pgas_m,
                                   "dens":self.pgas_dens,
                                   "met": 12+np.log10(self.pgas_met)})
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

    def AIC_test(self,dfin,breakpoint):
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

        med_data = medfilt(dfin['rad'], kernel_size=21)
        popt,pcov = curve_fit(UTILITY.linear_fit, dfin['rad'],dfin['met'])

        y1 = UTILITY.linear_fit(dfin['rad'],*popt)
        
        x0 = np.array([min(dfin['rad']), breakpoint, max(dfin['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(dfin['rad'], dfin['met'])
        my_pwlf.fit_with_breaks(x0)
        
        xHat = np.linspace(min(dfin['rad']), max(dfin['rad']), num=10000)
        yHat = my_pwlf.predict(xHat)
    

        RSS_linear = np.sum((dfin['met']-(y1))**2)
        RSS_broken = np.sum((dfin['met']-(yHat))**2)
        print(RSS_broken)
        print(RSS_linear)
        AIC_L = 4+(len(y1)*np.log(RSS_linear))
        AIC_B = 8+(len(yHat)*np.log(RSS_broken))

        print("AIC VALUES: \n ")
        print("Linear {}".format(AIC_L))
        print("Broken {}".format(AIC_B))
        

        if AIC_L>AIC_B:
            val = 0
        elif AIC_B>AIC_L:
            val=1
        return val

    def fit_linear(self,dfin, pc):
        '''
        Pseudocode
        INPUTS:
        - dfin -> subhalo dataframe
        - pc -> annuli pc (0.1->1)
        Calculate linear fit (f(x) = a*x+b)
        Take input of dataframe and fit linear trendline using scipy curve_fit 
        (popt,pcov) (popt[0] = gradient) , (popt[1]=intercept)

        return -> plot png file saved to directory with ID, simulation and snapshot as filename/title
        '''
        annuli = pc
        dfin = dfin[dfin['rad']<pc]     
        popt,pcov = curve_fit(UTILITY.linear_fit, dfin['rad'],dfin['met'])
        med_data = medfilt(dfin['met'],kernel_size = 21)

        plt.figure(figsize=(20,12))
        plt.plot(dfin['rad'], med_data, 'r-')
        plt.plot(dfin['rad'], UTILITY.linear_fit(dfin['rad'],*popt))
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10} (\frac{O}{H})$")
        filename = 'Visuals/linearfit/sub_{}_met_snap={}'.format(self.subID,self.snapID)
        plt.savefig(filename)
        plt.close()
    
    def fit_quad(self,dfin,pc):
        '''
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
        annuli = pc
        dfin.sort_values(by='rad',inplace = True)
        med_data1 = medfilt(dfin['met'], kernel_size=21)
        x0 = np.array([min(dfin['rad']), breakpoint, max(dfin['rad'])])
        break1 = list(dfin[dfin['rad']<breakpoint])
        break2 = list(dfin[dfin['rad']>=breakpoint])
        my_pwlf = pwlf.PiecewiseLinFit(dfin['rad'], dfin['met'])
        my_pwlf.fit_with_breaks(x0)
        
        xHat = np.linspace(min(dfin['rad']), max(dfin['rad']), num=10000)
        yHat = my_pwlf.predict(xHat)
        '''
        popt1,pcov1 = curve_fit(UTILITY.linear_fit, break1['rad'],break1['met'])[0]
        popt2,pcov2 = curve_fit(UTILITY.linear_fit, break2['rad'],break2['met'])[0]
        
        p1 = (popt1*break1)+pcov1
        p1=np.append(p1,((popt2*break2)+pcov2))
        
        '''
        
        plt.figure(figsize=(20,12))
        plt.plot(dfin['rad'], med_data1, 'b--')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        filename = 'Visuals/breakfit/33/sub_{}_snapshot_{}.png'.format(self.subID, self.snapID)
        plt.savefig(filename)
        plt.close()
    
    def savgol_smooth(self,dfin,pc,scatter_q):
        '''
        Pseudocode

        Takes input of dataframe with 'rad' and 'met' keywords, which represent radius and metlalicity 
        Calculates a Savitzky-Golay filter -> used to smooth metallicity data for clearer representation of the gradient 

        WARN - values calculated by the Savitzky-Golay filter cannot be used to calculate scientific results 

        return -> plot png file saved to directory with ID, simID and snapID as filename/title

        Optional Extras -> could include conditions to plot additional information onto graph such as linear/broken/quadratic fits? 
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
        '''
        Pseudocode 

        Takes input of dataframe and applies bootstrapping depending on type (either linear or quadratic)

        samples dataframe keeping only fraction of values (random) each time -> 

        compare means, range 

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
 
    def gas_visualisation(self, dfin, decp,pc):
        '''
        Psuedocode
        Inputs -> dataframe containing gas density data (aligned to z axis)

        flattens gas density data by decp (0 = 1kpc box, 1 .1kpc , 2 = .01kpc etc)
        
        plot visual by scatter of density values -> see other codes for plot syntax
        
        '''
        df = dfin.round(decp)
        annuli = pc
        df= df[df['rad']<annuli]
        df = df.groupby(['x','y'])['dens'].sum().reset_index()
        plt.figure(figsize=(20,12), dpi=500)
        plt.style.use('dark_background')
        plt.scatter(df['x'],-df['y'],c=(np.log10(df['m'])),cmap='inferno', vmin=(min(np.log10(df['m']))),vmax =(max(np.log10(df['m']))))
        plt.xlabel('$\Delta x$ [kpc/h]')
        plt.ylabel('$\Delta y$ [kpc/h]')
        plt.colorbar(label='log10(Gas Mass)')
        plt.title('Gas Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
        filename = 'Visuals/visuals/Mgass_{}_sub_{}.png'.format(self.simID, self.subID)
        plt.savefig(filename)
        plt.close()
    
    def gen_slope(self,dfin):
        df = dfin
        popt,pcov = curve_fit(UTILITY.linear_fit, df['rad'],df['met'])
        met = self.tot_met
        mass = self.m_tot
        sfr = self.totsfr
        return (popt[0],met,mass,sfr)


'''
sub = galaxy("TNG50-1",99,117260)
sub.galcen()
sub.ang_mom_align('gas')
sub.rad_transform()
dfg = sub.df_gen('gas','comb')
dfg2 = sub.rad_norm(dfg,10)
dfg2 = sub.z_filter(dfg2)
sub.AIC_test(dfg2,3)
#sub.broken_fit(dfg2,3,10)
#sub.savgol_smooth(dfg2,10,'Y')
'''

dfin = pd.read_csv("csv/tng99subhalos.csv")
valid_id= list(dfin['id'])
def slopeplot_dataget(i):
    try:
        sub = galaxy("TNG50-1",99,i)
        sub.galcen()
        sub.ang_mom_align('gas')
        sub.rad_transform()
        dfg = sub.df_gen('gas','comb')
        dfg2 = sub.rad_norm(dfg,10)
        dfg2 = sub.z_filter(dfg2)
        AICval = sub.AIC_test(dfg2,0.4)
        slope,met,mass,sfr = sub.gen_slope(dfg2)
        id = i
        print("subhalo {} calculated".format(i))
        return (slope,met,mass,id,sfr,AICval)
    except OSError:
        return print('OSerror')
    '''
    except ValueError:
        return print("value-error")
    except KeyError:
        return print("keyerror")
    except OSError:
        return print('OSerror')
    except TypeError:
        return print('TypeError')
    '''
#0 = linear>broken
#1 = broken>linear
xval = np.linspace(0,13,100)
def line(m,x,b):
    y = 10**((m*x)+b)
    return y 

returns = Parallel(n_jobs=25)(delayed(slopeplot_dataget)(i) for i in valid_id)
df2=pd.DataFrame(returns,columns=['slope','met','mass','id','sfr','AICval'])
df2.to_csv("csv/tng99slopes.csv")

end = time.time()
print('runtime = {}s'.format(end-start))
