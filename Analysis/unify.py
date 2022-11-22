import time  # runtime calculation import numpy as np #data handling
import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests  # obtain data from API server
import pwlf
from joblib import Parallel, delayed
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter
import BCUTILS
from BCUTILS import UTILITY
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()


class subhalo:
    def __init__(self,simID,snapID,subID):
        '''
        Subhalo Class -> repeat creation of subhalo object and analyis tools to determine properties, fit gradients and statistical analysis 
         
        
        '''
        basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'
        baseurl = 'http://www.tng-project.org/api/'+str(simID)+'/snapshots/'+str(snapID)
        
        hubble = 0.7
        #check snapID value to obtain correct redshift and scale factor 
        if snapID is 99:
            redshift = 0.0
        elif snapID is 67:
            redshift =0.503047523244883
        elif snapID is 33:
            redshift = 2.00202813925285
        else:
            data = BCUTILS.get(baseurl)
            redshift = data['redshift']
        self.redshift = redshift
        scalefac = 1./(1.+redshift) #calculate scale factor
        ptNumGas = il.snapshot.partTypeNum('gas') #obtain part index value to properly query TNG data
        ptNumStars = il.snapshot.partTypeNum('stars')
        #
        # pull global data 
        #
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
        
        #
        # pull gas particle level data 
        #
        
        gas = il.snapshot.loadSubhalo(basePath, snapID, subID, 'gas', fields=['Coordinates', 'Masses','Density','Velocities', 'StarFormationRate','GFM_Metallicity'])
        # dimensions and units (see https://www.tng-project.org/data/docs/specifications/#parttype0):
        # Coordinates (N,3) ckpc/h   where ckps stands for co-moving kpc
        # Masses      (N)   10**10 Msun/h
        # Velocities  (N,3) km sqrt(scalefac)        # We convert these to pkpc (proper kpc), Msun and km/s, respectively
        crit_dist = 5 * self.Rhalf #30. # proper kpc
        self.crit_dist = crit_dist
        hcoldgas  = np.where( (gas['StarFormationRate'] > 0.0) & (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        #hcoldgas  = (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2)
        self.pgas_coo   = gas['Coordinates'][hcoldgas]/hubble / (1. + redshift)
        self.pgas_m     = gas['Masses'][hcoldgas] * 10**10 / hubble
        self.pgas_vel   = (gas['Velocities'][hcoldgas] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.conv_kms2kpcyr = (3.1558 / 3.08568) * 10**(-9)
        self.pgas_vel   = self.pgas_vel * self.conv_kms2kpcyr    #Convert to kpc/yr
        self.pgas_sfr   = gas['StarFormationRate'][hcoldgas]
        self.pgas_met   =gas['GFM_Metallicity'][hcoldgas]
        self.pgas_dens = gas['Density'][hcoldgas]
        
        #
        #pull stellar particle levl data 
        #
        stars = il.snapshot.loadSubhalo(basePath, snapID, subID, 'stars', fields=['Coordinates', 'Masses', 'Velocities','GFM_Metallicity' ])
        hstar = np.where( (np.sum((stars['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        self.pstar_coo   = stars['Coordinates'][hstar]/hubble / (1. + redshift)
        self.pstar_m     = stars['Masses'][hstar] * 10**10 / hubble
        self.pstar_vel   = (stars['Velocities'][hstar] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.pstar_vel   = self.pstar_vel * self.conv_kms2kpcyr
        self.pstar_met = stars['GFM_Metallicity'][hstar]
        
    def galcen(self):
        '''
        simple function to centre subhalo on centre of mass co-ordinates (normalisation of sorts)
        '''
        self.pgas_coo -= self.centre[None,:]
        self.pstar_coo -= self.centre[None,:]
        
    def ang_mom_align(self, type):
        '''
        Function to read particle level data -> using velocities, masses and positions calculate angular velocities and create transformation matrix to 
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
                                   "met":(self.pgas_met),
                                   "met2":(self.pgas_met)})
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
    
    def df_process(self,dfin,scale):
        '''
        Alternative to running rad_norm and z_filter if both being used 
        '''
        df = dfin
        df.rad = scale*((df.rad-df.rad.min())/(df.rad.max()-df.rad.min()))
        scaleheight = 0.1* self.stellarphotometricsrad
        dfout = df[df['z']<scaleheight]
        self.df = dfout
        return dfout
    
    def slopegen(self,breakpoint):
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
        df = self.df
        popt,pcov = curve_fit(UTILITY.linear_fit, df['rad'],df['met'])

        y1 = UTILITY.linear_fit(df['rad'],*popt)
        AIC = 1
        x0 = np.array([min(df['rad']), breakpoint, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], df['met'],disp_res = True)
        brokenslope = my_pwlf.slopes[0]
        my_pwlf.fit_with_breaks(x0)
        
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
        
        if linear>broken:
            return brokenslope
        elif broken>linear:
            return popt[0]
        else:
            return popt[0]

#
#Moving to more accessory functions -> graph plotting with fitting and gas density visualisation
#
    
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
        r'''
        Pseudocode 
        
        Inputs:
        
        dfin: pandas dataframe
            Pandas dataframe containing subhalo information -> function df_gen in this subhalo class provides correct keys
            Must contain ['met','rad']             
        breakpoint: int,
            Radial value for the breakpoint (the point at which break in linear fit is placed)

        Process                             
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
        plt.savefig('test2.png')
        plt.close()
    
    def gas_visualisation(self, dfin, decp):
        '''
        Psuedocode
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
        
ids = UTILITY.get_ids(99)
print(ids)