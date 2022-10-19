
import logging # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}

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

'''
def geturl(simID, snapID):
    baseurl = 'http://www.tng-project.org/api/'+str(simID)+'/snapshots/'+str(snapID)
    return baseurl
'''

class galaxy:
    def __init__(self,simID,snapID,subID):
        self.simID = simID
        self.snapID = snapID
        self.subID = subID 
        basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output'
        baseurl = 'http://www.tng-project.org/api/'+str(simID)+'/snapshots/'+str(snapID)
        
        hubble =0.7
        
        snap_util =get(baseurl) #use api to easily read snapshot information
        redshift = snap_util['redshift']
        self.redshift =redshift 
        
        scalefac = 1./(1.+redshift)
        
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
        # [units: proper kpc] (quantified in 3D)
        
        # Load all the relevant particle level info
        gas = il.snapshot.loadSubhalo(basePath, snapID, subID, 'gas', fields=['Coordinates', 'Masses', 'Velocities', 'StarFormationRate','GFM_Metallicity'])
        # dimensions and units (see https://www.tng-project.org/data/docs/specifications/#parttype0):
        # Coordinates (N,3) ckpc/h   where ckps stands for co-moving kpc
        # Masses      (N)   10**10 Msun/h
        # Velocities  (N,3) km sqrt(scalefac)        # We convert these to pkpc (proper kpc), Msun and km/s, respectively
        crit_dist = 5 * self.Rhalf #30. # proper kpc
        #np.where( (gas['StarFormationRate'] > 0.) & 
        hcoldgas  = (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2)
        self.pgas_coo   = gas['Coordinates'][hcoldgas]/hubble / (1. + redshift)
        self.pgas_m     = gas['Masses'][hcoldgas] * 10**10 / hubble
        self.pgas_vel   = (gas['Velocities'][hcoldgas] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.conv_kms2kpcyr = (3.1558 / 3.08568) * 10**(-9)
        self.pgas_vel   = self.pgas_vel * self.conv_kms2kpcyr    #Convert to kpc/yr
        self.pgas_sfr   = gas['StarFormationRate'][hcoldgas]
        self.pgas_met   =gas['GFM_Metallicity']
        
        # Load all stellar particle data
        stars = il.snapshot.loadSubhalo(basePath, snapID, subID, 'stars', fields=['Coordinates', 'Masses', 'Velocities'])
        hstar = np.where( (np.sum((stars['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        self.pstar_coo   = stars['Coordinates'][hstar]/hubble / (1. + redshift)
        self.pstar_m     = stars['Masses'][hstar] * 10**10 / hubble
        self.pstar_vel   = (stars['Velocities'][hstar] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.pstar_vel   = self.pstar_vel * self.conv_kms2kpcyr
        
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
        
    def visualise_cutout_TNG50(id, type, lim):


        sub_prog_url = "http://www.tng-project.org/api/TNG50-2/snapshots/99/subhalos/"+str(id)+"/"
        sub_prog = get(sub_prog_url)
        
        cutout_request = {'gas':'Coordinates,Masses,GFM_Metallicity'}
        cutout = get(sub_prog_url+"cutout.hdf5", cutout_request)


        with h5py.File(cutout,'r') as f:
            x = f['PartType0']['Coordinates'][:,0] - sub_prog['pos_x']
            y = f['PartType0']['Coordinates'][:,1] - sub_prog['pos_y']
            if (type=='gas'):
                dens =np.exp(f['PartType0']['Masses'][:])**4
            elif(type =='met'): 
                dens =-np.log10(f['PartType0']['GFM_Metallicity'][:])**3

        plt.figure()
        plt.hist2d(x,y,weights=dens,bins=[1500,1000], cmap = 'afmhot', vmin = min(dens), vmax = (max(dens)))
        plt.xlabel('$\Delta x$ [ckpc/h]')
        plt.ylabel('$\Delta y$ [ckpc/h]')
        plt.xlim(-10,10)
        plt.ylim(-7,4)
        plt.savefig('hist_{}_{}.png'.format(type,id))
        plt.close()

        return print("graph plotted for subhalo{}".format(id))
    def visualise_cutout_TNG100(id, type, lim):


        sub_prog_url = "http://www.tng-project.org/api/TNG100-1/snapshots/99/subhalos/"+str(id)+"/"
        sub_prog = get(sub_prog_url)
        
        cutout_request = {'gas':'Coordinates,Masses,GFM_Metallicity'}
        cutout = get(sub_prog_url+"cutout.hdf5", cutout_request)


        with h5py.File(cutout,'r') as f:
            x = f['PartType0']['Coordinates'][:,0] - sub_prog['pos_x']
            y = f['PartType0']['Coordinates'][:,1] - sub_prog['pos_y']
            if (type=='gas'):
                dens =np.exp(f['PartType0']['Masses'][:])**4
            elif(type =='met'): 
                dens =-np.log10(f['PartType0']['GFM_Metallicity'][:])**3

        plt.figure()
        plt.hist2d(x,y,weights=dens,bins=[1500,1000], cmap = 'afmhot', vmin = min(dens), vmax = (max(dens)))
        plt.xlabel('$\Delta x$ [ckpc/h]')
        plt.ylabel('$\Delta y$ [ckpc/h]')
        plt.xlim(-10,10)
        plt.ylim(-7,4)
        plt.savefig('hist_{}_{}.png'.format(type,id))
        plt.close()

        return print("graph plotted for subhalo{}".format(id))
    zval =np.empty(shape=(2,2))
    
    #def sub_flatten(self):
    #    zlen = len(self.pgas_coo[:,2])
    #    xlist = self.pgas_coo[:,0]
    #    ylist = self.pgas_coo[:,1]
    #    for i,j in range(zlen):
    #        = np.sum(self.pgas_coo[i][j][])
    


sub0 = galaxy('TNG100-1',70,0)
sub0.galcen()
sub0.ang_mom_align('gas')
'''
print(sub0.pgas_coo)
print(np.ndim(sub0.pgas_m))
'''
x= sub0.pgas_coo[:,0]
y= sub0.pgas_coo[:,1]
dens = (sub0.pgas_m)
plt.figure()
plt.scatter(x,y,c=dens,cmap='inferno')
plt.savefig('pls.png')
plt.close()
df = pd.DataFrame({"x": sub0.pgas_coo[:,0],"y": sub0.pgas_coo[:,1],"z": sub0.pgas_coo[:,2],"m": sub0.pgas_m})
df2=df.round()
df2.sort_values(by='x',inplace=True)
df2.to_csv('filet.csv')
