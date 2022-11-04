#%%
from itertools import groupby
import logging
from random import random # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il
from scipy.signal import medfilt
from joblib import Parallel, delayed
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
from sklearn.preprocessing import MinMaxScaler
#basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output'
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

class galaxy:
    def __init__(self,simID,snapID,subID):
        
        #object creation requires only 3 input parameters for galaxy selection 
        self.simID = simID #simulation used (TNG100/TNG50)
        self.snapID = snapID #snapshot of study - can be either number of z= (z= required in parenthesis )
        self.subID = subID #subhalo ID - could be iterated 
        #basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output' #path to simulation data on Noether
        basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/'
        baseurl = 'http://www.tng-project.org/api/'+str(simID)+'/snapshots/'+str(snapID) # API url for simulation/snapshot information
        
        hubble =0.7
        
        snap_util =get(baseurl) #use api to easily read snapshot information
        redshift = snap_util['redshift']
        self.redshift =redshift  # store redshift value as attribute
        
        scalefac = 1./(1.+redshift) #calculate scale factor
        
        #
        # Read Subhalo level info 
        #

        ptNumGas = il.snapshot.partTypeNum('gas') #determine index designation for each particle type
        ptNumStars = il.snapshot.partTypeNum('stars')
        #pull all data for specific subhalo 
        all_fields= il.groupcat.loadSingle(basePath, snapID, subhaloID = subID)
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
        hcoldgas  = np.where( (gas['StarFormationRate'] > 0.) & (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
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
    def radial_calc(self):
        self.pgas_rad_len   = np.sqrt(self.pgas_coo[:,0]**2+self.pgas_coo[:,1]**2)
        self.pstar_rad_len  = np.sqrt(self.pstar_coo[:,0]**2+self.pstar_coo[:,1]**2)
        #print(self.pgas_rad_len)
        #pgas_rad=[]
        #numvals = np.linspace(0,len(self.pgas_coo[:,0]),1)
        '''
        for i in numvals:
            pgas_rad.append(np.sqrt(self.pgas_coo[:,0][i]**2+self.pgas_coo[:,1][i]**2))
        rmax = max(pgas_rad)
        rmin = min(pgas_rad)
        rnorm=[]
        for i in pgas_rad:
            rnorm.append((100*(pgas_rad[i]-rmin)/(rmax-rmin)))
        self.rnorm = rnorm
        '''



    def dataframegen(self,type):
        if(type=='gas'):
            df=pd.DataFrame({"x": self.pgas_coo[:,0], "y":self.pgas_coo[:,1],"z":self.pgas_coo[:,2], "rad":abs(self.pgas_rad_len),"m":self.pgas_dens,"met": self.pgas_met})
            #print(df['rad'])
        elif(type=='star'):
            df=pd.DataFrame({"x": self.pstar_coo[:,0], "y":self.pstar_coo[:,1],"z":self.pstar_coo[:,2], "rad":abs(self.pstar_rad_len),"m":self.pstar_m,"met":self.pstar_met})
        self.df = df
        return df

    #'''
    def height_filter(self,dfin):
        z_max = 0.1*self.stellarphotometricsrad
        df_filtered = dfin[dfin['z']<z_max]
        self.filt_df = df_filtered
        return df_filtered
    '''
    def rad_normalise(self,df):
        rad_arr = df['rad'].to_numpy()
        rmax = max(rad_arr)
        rmin = min(rad_arr)
        numvals = np.arange(0,len(rad_arr),1)
        rad2 = []
        for j in numvals:
            rad2.append(100*(rad_arr[j]-rmin)/(rmax-rmin))
        rad2df = pd.DataFrame({'rad':rad2})
        df['rad'] = rad2df['rad'].values
        return df
    '''

        


class visualisation:
    def __init__(self, df_g, df_s, subID, snapID, simID,crit_dist):
        self.df_g = df_g #gas dataframe
        self.df_s = df_s #stars dataframe
        self.subID = subID
        self.snapID = snapID
        self.simID = simID
        self.crit_dist = crit_dist
    def visual(self,type, quant, decp, annuli_pc):
        if (type=='gas'):
            df = self.df_g
            if(quant=='mass'):
                df_valid = df.round(decp)
                annul1= annuli_pc*self.crit_dist
                df_valid = df[df['rad']<annul1]
                df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
                plt.figure(figsize=(20,12), dpi=500)
                plt.style.use('dark_background')
                plt.scatter(-df_valid['x'],-df_valid['y'],c=(np.log10(df_valid['m'])),cmap='inferno', vmin=(min(np.log10(df_valid['m']))),vmax =(max(np.log10(df_valid['m']))))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='log10(Gas Mass)')
                plt.title('Gas Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
                filename = 'suppng/Mgass_{}_sub_{}.png'.format(self.simID, self.subID)
                plt.savefig(filename)
                plt.close()
            elif(quant=='metallicity'):
                df_valid = df.round(decp)
                annul1= annuli_pc*self.crit_dist
                df_valid = df[df['rad']<annul1]
                df_valid = df_valid.groupby(['x','y'])['met'].sum().reset_index()
                plt.figure(figsize=(20,12), dpi=500)
                plt.style.use('dark_background')
                plt.scatter(-df_valid['x'],-df_valid['y'],c=(np.log10(df_valid['m'])),cmap='inferno', vmin=(min(np.log10(df_valid['m']))),vmax =(max(np.log10(df_valid['m']))))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='log10(Gas Metallicity)')
                plt.title('Metallicity Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
                filename = 'temppng/met_gas_{}_sub_{}.png'.format(self.simID, self.subID)
                plt.savefig(filename)
                plt.close()

        elif(type=='stars'):
            df = self.df_s
            if (quant=='mass'):
                df_valid = df.round(decp)
                annul1= annuli_pc*self.crit_dist
                df_valid = df[df['rad']<annul1]
                df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
                plt.figure(figsize=(21,15))
                plt.style.use('dark_background')
                plt.scatter(-df_valid['x'],-df_valid['y'],c=(np.log10(df_valid['m'])),cmap='inferno', vmin=np.log10(min(df_valid['m'])), vmax = np.log10(max(df_valid['m'])))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='log10(Stellar Mass)')
                plt.title('Gas Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
                filename = 'temppng/Mstar_{}_sub_{}.png'.format(self.simID, self.subID)
                plt.savefig(filename)
                plt.close()
            
            elif(quant=='metallicity'):
                df_valid = df.round(decp)
                annul1= annuli_pc*self.crit_dist
                df_valid = df[df['rad']<annul1]
                df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
                plt.figure(figsize=(21,15))
                plt.style.use('dark_background')
                plt.scatter(-df_valid['x'],-df_valid['y'],c=(df_valid['met']),cmap='inferno', vmin=min(df_valid['m']), vmax = max(df_valid['m']))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='Metallicity)')
                plt.title('Metallicity Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
                filename = 'temppng/met_star_{}_sub_{}.png'.format(self.simID, self.subID)
                plt.savefig(filename)
                plt.close()

    def metgrad (self,type,decp, annuli_pc):
        if (type=='gas'):
            df = self.df_g
            annul1= annuli_pc*self.crit_dist
            df_valid = df[abs(df['rad'])<annul1]
            df_valid.rad = 100*((df_valid.rad-df_valid.rad.mean())/(df_valid.rad.max()-df_valid.rad.min()))
            df_valid.sort_values(by='rad',inplace=True)
            df_valid = df_valid[abs(df_valid['met'])>0]
            medfit = medfilt(df_valid['met'],kernel_size=21)
            plt.figure(figsize=(21,15))
            plt.plot(df_valid['rad'], (12+np.log10(medfit)), 'r--')
            #plt.scatter(df_valid['rad'], (12+np.log10(df_valid['met'])), c=df_valid['m'], cmap = 'viridis')
            plt.xlabel('Radial Distance [kpc/h]')
            plt.ylabel('12+log10(O/H) [kpc/h]')
            #plt.colorbar(label='Gass mass')
            plt.title('Gas Metallicity Gradient for SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
            filename = 'metgradGAS_{}_sub_{}.png'.format(self.simID, self.subID)
            plt.savefig(filename)
            plt.close()

        elif (type=='stars'):
            df = self.df_s
            annul1= annuli_pc*self.crit_dist
            df_valid = df[df['rad']<annul1]
            df_valid=df_valid.round()
            df_valid.groupby(['rad'])['met'].mean().reset_index()
            df_valid.sort_values(by='rad', inplace=True)
            plt.figure(figsize=(21,15))
            plt.scatter(df_valid['rad'], (12+np.log10(df_valid['met'])), c=df_valid['m'], cmap = 'viridis')
            plt.xlabel('Radial Distance [kpc/h]')
            plt.ylabel('12+log10(O/H) [kpc/h]')
            plt.colorbar(label='Gass mass')
            plt.title('Stellar Metallicity gradient for SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
            #temppng/
            filename = 'metgradSTAR_{}_sub_{}.png'.format(self.simID, self.subID)
            plt.savefig(filename)
            plt.close()
    def  plot3d(self):
        df = self.df_g
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        ax.scatter(df['x'], df['y'], df['z'], cmap = 'inferno')
        plt.savefig("3d.png")
        plt.show()
        plt.close()





galaxy_df = pd.read_csv("test1.csv")
valid_galaxies = galaxy_df[galaxy_df['mass']<9.5]

'''
sub1 = galaxy("TNG50-1", 99, 307614)
sub1.galcen()
sub1.ang_mom_align('gas')
sub1.radial_calc()
dfg = sub1.dataframegen('gas')
dfs = sub1.dataframegen('star')
dfg2=sub1.height_filter(dfg)
dfs2=sub1.height_filter(dfs)
sub1plot = visualisation(dfg2, dfs2,sub1.subID, sub1.snapID, sub1.simID, sub1.crit_dist)
sub1plot.metgrad('gas',1,1)
'''
#'''
def met_plot(i):
    sub1 = galaxy("TNG50-1", 99, i)
    sub1.galcen()
    sub1.ang_mom_align('gas')
    sub1.radial_calc()
    dfg = sub1.dataframegen('gas')
    dfs = sub1.dataframegen('star')
    dfg2=sub1.height_filter(dfg)
    dfs2=sub1.height_filter(dfs)
    #dfg2 = sub1.rad_normalise(dfg2)
    sub1plot = visualisation(dfg2, dfs2,sub1.subID, sub1.snapID, sub1.simID, sub1.crit_dist)
    sub1plot.plot3d()
    return print(".") 
met_plot(3052)
#returns = Parallel(n_jobs=15)(delayed(met_plot)(i) for i in valid_galaxies['id'])
#'''

end = time.time()
print("runtime {}".format(end-start))

# %%
