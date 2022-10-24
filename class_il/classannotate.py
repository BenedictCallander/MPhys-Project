#import classes 
import logging # http logging for debugging purpouses
from random import random #for dataset thinning in debug 
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd #data storage and manipulation
import numpy as np #mathematical functions
import matplotlib.pyplot as plt #plotting libary 
import illustris_python as il #illustris file functions


headers = {"api-key":"YOUR-API-KEY"}
print("start time: {}".format(time.time())) #print start time
start =time.time() #runtime counter start

#basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output'
#utility function to read data from Illustris API 

def get(path, params = None):

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

#galaxy class - to create identically formatted objects and 
#contain utility functions 
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

        ptNumGas = il.snapshot.partTypeNum('gas') #determine index designation for each particle type (gas=0 stars = 4)
        ptNumStars = il.snapshot.partTypeNum('stars')
        #pull all data for specific subhalo 
        all_fields= il.groupcat.loadSingle(basePath, snapID, subhaloID = subID)
        self.lMgas  = np.log10( all_fields['SubhaloMassInRadType'][ptNumGas]/hubble ) + 10.
        self.lMstar = np.log10( all_fields['SubhaloMassInRadType'][ptNumStars]/hubble ) + 10.
        # Coordinate of particle with minimum binding energy (converted from ckpc/h to kpc)
        self.centre = all_fields['SubhaloPos']/hubble / (1. + redshift)  # 3-element array [units: proper kpc]
        # Adopt the 3D half-stellar-mass radius
        self.Rhalf  = all_fields['SubhaloHalfmassRadType'][ptNumStars]/hubble / (1. + redshift)  
        # [units: proper kpc] (quantified in 3D)
        
        # Load all the relevant particle level info
        gas = il.snapshot.loadSubhalo(basePath, snapID, subID, 'gas', fields=['Coordinates', 'Masses', 'Velocities', 'StarFormationRate','GFM_Metallicity'])
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

    def radial_coo(self):
        self.pgas_rad_len   = np.sqrt(self.pgas_coo[:,0]**2+self.pgas_coo[:,1]**2)
        self.pstar_rad_len  = np.sqrt(self.pstar_coo[:,0]**2+self.pstar_coo[:,1]**2)
    def xyzprop_df(self,type):
        if (type=='met'):
            df=pd.DataFrame({"x": self.pstar_coo[:,0], "y":self.pstar_coo[:,1],"z":self.pstar_coo[:,2], "rad":self.pstar_rad_len,"m": 10e9*self.pstar_met})
        elif (type=='mass'):
            df=pd.DataFrame({"x": self.pstar_coo[:,0], "y":self.pstar_coo[:,1],"z":self.pstar_coo[:,2], "rad":self.pstar_rad_len,"m": self.pstar_m})
        
        self.df = df
        return df # in case dataframe object wanted as return 
    def gas_plot(self,annuli_pc):
        
        annul1= annuli_pc*self.crit_dist
        df_valid = self.df[self.df['rad']<annul1]
        df_valid =df_valid.round(3)
        df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
        plt.figure(figsize=(21,15))
        plt.style.use('dark_background')
            #plt.hist2d(data2['x'],data2['y'], weights=data2['m'],bins=[500,500],cmap = 'inferno',vmin=(min(data2['m'])),vmax = max(data2['m']))
        plt.scatter(-df_valid['x'],-df_valid['y'],c=(np.log10(df_valid['m'])),cmap='inferno', vmin=min(np.log10(df_valid['m'])), vmax = max(np.log10(df_valid['m'])))
            #plt.plot(critical_radius,zeros,'g+', markersize=10,label='critical radius')
        plt.xlabel('$\Delta x$ [kpc/h]')
        plt.ylabel('$\Delta y$ [kpc/h]')
            #plt.xlim(-100,100)
            #plt.ylim(-100,100)
        plt.colorbar(label='log10(Gas Density)')
        plt.title('Gas Density of SubID {}: TNG100-1 snapshot 70'.format(self.subID))
            #plt.legend(loc='upper right')
        filename = 'temppng/TNG50_sub_{}.png'.format(sub1.subID)
        plt.savefig(filename)
        plt.close()
    def met_grad_gas(self,annuli_pc):
        df=pd.DataFrame({"x": self.pgas_coo[:,0], "y":self.pgas_coo[:,1],"z":self.pgas_coo[:,2], "rad":self.pgas_rad_len,"metallicity": self.pgas_met})
        annul1= annuli_pc*self.crit_dist
        met_df = df[df['rad']<annul1]
        
        plt.figure(figsize=(21,15))
        plt.plot(met_df['rad'], 12+np.log10(met_df['metallicity']), 'b.')
        plt.xlabel('X')
        plt.ylabel('metallicity')
        filename = 'temppng/TNG50_metgrad_sub_{}.png'.format(sub1.subID)
        plt.savefig(filename)
        plt.close()
        return print('.')
    def met_grad_star(self,annuli_pc):
        df =pd.DataFrame({"rad": self.pstar_rad_len,"met": self.pstar_met })
        #df= df.sample(frac=0.01,random_state=1)
        annul1= annuli_pc*self.crit_dist
        met_df = df[df['rad']<annul1]
        plt.figure(figsize=(21,15))
        plt.scatter(met_df['rad'], 12+np.log10(met_df['met']), 'b.')
        plt.xlabel('Radial Position')
        plt.ylabel('metallicity')
        filename = 'temppng/TNG50_metgrad__stellar_forsub_{}.png'.format(sub1.subID)
        plt.savefig(filename)
        plt.close()
        return print('.')

'''
unrequired function for finding most massive star forming subhalos in simulation 
to find for other than TNG 100-1 snap 70:
alter:
-massive_url TNG value 
- snapshot value 

massive_url = "http://www.tng-project.org/api/TNG100-1/snapshots/70/subhalos/?order_by=-mass&sfr_gt=0.0/"

for i in range (20):
    valid_subs = get(massive_url)
    massive_ids = [ valid_subs['results'][i]['id'] for i in range(20)]
print(massive_ids)
#massive_ids = [0, 7516, 21013, 15129, 31129, 39628, 47416, 26558, 44002, 34668, 57620, 51083, 54570, 63544, 69982, 60421, 85032, 87479, 88730, 80680]
'''
massive_list=[0, 63864, 96762, 117250, 143880, 184931, 198182, 208811, 220595, 229933, 253861, 167392, 242788, 264883, 282779, 275545, 294866, 289385, 313692, 300903]
#old debug print statements - returns length of pgas co-ordinates and metallicity (equiv to number of cells in subhalo)
#print(len(sub1.pgas_coo))
#print(len(sub1.pgas_met))

#for i in massive_ids:
#iterative loop for a range of subhalo IDS  
sub1= galaxy('TNG50-1',99,63864) #create galaxy object using galaxy class - and arguments for __init__ 
sub1.galcen() #centre subhalo 
sub1.ang_mom_align('gas') #align subhalo so is parallel with z-axis - can specify either gas or stars 
sub1.radial_coo() #calculate radial co-ordinates
sub1.xyzprop_df('mass') #generate self.df - specifiers are mass or met for gass mass or metallicity 
sub1.gas_plot(0.25) #gas visualisation, argument is % of critical radius to be used for annulus 


end = time.time()
print("programme end time :{}".format(time.time()))
print('runtime = {} seconds'.format(end-start))