import logging # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il
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


class galaxy(object):
    def __init__(self, sim_name='TNG100-1', snapNum=70, subhaloID=0):
        self.sim_name  = sim_name
        self.snapNum   = snapNum

        basePath = '/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/' + sim_name + '/output'
        baseUrl  = 'http://www.tng-project.org/api/' + sim_name '/snapshots/1'
        hubble   = 0.7

        # Determine the redshift and scale factor of the Universe corresponding to snapNum (i.e. epoch of observation)
        snap_list = get(baseUrl)
        redshift = snap_list[snapNum]['redshift']
        scalefac = 1. / (1. + redshift)
        self.redshift = redshift
        
        # Load all the relevant subhalo level info
        ptNumGas  = il.snapshot.partTypeNum('gas') # 4
        ptNumStars  = il.snapshot.partTypeNum('stars') # 4
        all_fields  = il.groupcat.loadSingle(basePath, snapNum, subhaloID=subhaloID) 
        self.lMgas  = np.log10( all_fields['SubhaloMassInRadType'][ptNumGas]/hubble ) + 10.
        self.lMstar = np.log10( all_fields['SubhaloMassInRadType'][ptNumStars]/hubble ) + 10.
        # Coordinate of particle with minimum binding energy (converted from ckpc/h to kpc)
        self.centre = all_fields['SubhaloPos']/hubble / (1. + redshift)  # 3-element array [units: proper kpc]
        # Adopt the 3D half-stellar-mass radius
        self.Rhalf  = all_fields['SubhaloHalfmassRadType'][ptNumStars]/hubble / (1. + redshift)  # [units: proper kpc] (quantified in 3D)
        
        # Load all the relevant particle level info
        gas = il.snapshot.loadSubhalo(basePath, snapNum, subhaloID, 'gas', fields=['Coordinates', 'Masses', 'Velocities', 'StarFormationRate'])
        # dimensions and units (see https://www.tng-project.org/data/docs/specifications/#parttype0):
        # Coordinates (N,3) ckpc/h   where ckps stands for co-moving kpc
        # Masses      (N)   10**10 Msun/h
        # Velocities  (N,3) km sqrt(scalefac)        # We convert these to pkpc (proper kpc), Msun and km/s, respectively
        crit_dist = 5 * self.Rhalf #30. # proper kpc
        hcoldgas  = np.where( (gas['StarFormationRate'] > 0.) & (np.sum((gas['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        self.pgas_coo   = gas['Coordinates'][hcoldgas]/hubble / (1. + redshift)
        self.pgas_m     = gas['Masses'][hcoldgas] * 10**10 / hubble
        self.pgas_vel   = (gas['Velocities'][hcoldgas] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.conv_kms2kpcyr = (3.1558 / 3.08568) * 10**(-9)
        self.pgas_vel   = self.pgas_vel * self.conv_kms2kpcyr    #Convert to kpc/yr
        self.pgas_sfr   = gas['StarFormationRate'][hcoldgas]
        
        # Load all stellar particle data
        stars = il.snapshot.loadSubhalo(basePath, snapNum, subhaloID, 'stars', fields=['Coordinates', 'Masses', 'Velocities'])
        hstar = np.where( (np.sum((stars['Coordinates']/hubble / (1. + redshift) - self.centre[None,:])**2, axis=1) < crit_dist**2) )[0]
        self.pstar_coo   = stars['Coordinates'][hstar]/hubble / (1. + redshift)
        self.pstar_m     = stars['Masses'][hstar] * 10**10 / hubble
        self.pstar_vel   = (stars['Velocities'][hstar] * np.sqrt(scalefac)) - all_fields['SubhaloVel'][None,:]
        self.pstar_vel   = self.pstar_vel * self.conv_kms2kpcyr


    #Calculating the centre of mass of a galaxy    
    def calc_galcen(self):

        # Uncomment if determining mass-weighted centre of gas particles
        self.pgas_coo  -= self.centre[None,:]
        self.pstar_coo -= self.centre[None,:]
        
    def calc_ang_mom(self, angmom_type='stars'):
        #== Compute the angular momentum vector ==
        # different possible angmom_type: gas', 'stars'
        if (angmom_type == 'gas'):
            _coo = np.copy(self.pgas_coo)
            _vel = np.copy(self.pgas_vel)
            _m   = np.copy(self.pgas_m)
        elif (angmom_type == 'stars'):
            _coo = np.copy(self.pstar_coo)
            _vel = np.copy(self.pstar_vel)
            _m   = np.copy(self.pstar_m)

        # Calculate angular momentum based on particle type of choice (stars, gas or baryons)
        self.ang_mom_3D       = np.sum(_m[:,None]*np.cross(_coo,_vel),axis=0) # (3-element array specifying orientation of angular momentum vector)
        self.ang_mom = self.ang_mom_3D/ np.sum(_m) #angular momentum vector divided by total mass 

        #Aligning the z axis to the angular momentum vector
        j = self.ang_mom/np.linalg.norm(self.ang_mom)  #normalised specific angular momentum vector

        x = np.array([1.,2.,3.]) #np.random.randn(3) # random vector # np.array([1.,1.,1.]) #
        x -= x.dot(j) * j # make it orthogonal to j
        x /= np.linalg.norm(x) # normalize it

        y = np.cross(j, x) #third vector orthogonal to x and j

        A=(x,y,j) #transformation matrix


        #== Work on gas particles ==
        self.pgas_coo=np.dot(A,self.pgas_coo.T).T  #change coordinates
        self.pgas_vel=np.dot(A,self.pgas_vel.T).T 
        
        #-- Make checkplot of x-y, y-z, x-z views --
        '''
        fig  = plt.figure()
        gs   = fig.add_gridspec(2, 2)

        idx  = np.random.uniform(low=0, high=len(self.pgas_coo[:,0]), size=5000).astype(int)#Draw random particles
        
        ax = fig.add_subplot(gs[0,0]) # x-y
        #ax.scatter(_coo[idx,0], _coo[idx,1], c="r", s=1)
        ax.scatter(self.pgas_coo[idx,0], self.pgas_coo[idx,1], s=1)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim((-10,10))
        ax.set_ylim((-10,10))
        
        ax = fig.add_subplot(gs[0,1]) # y-z
        #ax.scatter(_coo[idx,1], _coo[idx,2], c="r", s=1)
        ax.scatter(self.pgas_coo[idx,1], self.pgas_coo[idx,2], s=1)
        ax.set_xlabel('y')
        ax.set_ylabel('z')
        ax.set_xlim((-10,10))
        ax.set_ylim((-10,10))
        
        ax = fig.add_subplot(gs[1,0]) # x-z
        #ax.scatter(_coo[idx,0], _coo[idx,2], c="r", s=1)
        ax.scatter(self.pgas_coo[idx,0], self.pgas_coo[idx,2], s=1)
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.set_xlim((-10,10))
        ax.set_ylim((-10,10))

        #plt.title('Galaxy in three different planes before and after the coordinate transformation')
        fig.tight_layout()
        fig.savefig('../plots/'+str(self.snapNum)+'_'+str(self.subhaloID)+'_3D.png')
        #plt.show()
        plt.close()
        '''
        #--
        
        #Changing coordinates to cylindrical polar
        phi= np.arctan(self.pgas_coo[:,1]/self.pgas_coo[:,0]) #phi component 
        phi[(self.pgas_coo[:,0] < 0)] += np.pi # if x is negative, need to add pi
        r  = np.sqrt(self.pgas_coo[:,0]**2 + self.pgas_coo[:,1]**2) #r component from cartesian to cylindrical
        v_r=(self.pgas_coo[:,0]*self.pgas_vel[:,0]+self.pgas_coo[:,1]*self.pgas_vel[:,1])/(np.sqrt(self.pgas_coo[:,0]**2 + self.pgas_coo[:,1]**2)) #radial velocity using the velocity and position in cartesian coordinates 
        v_phi=(self.pgas_coo[:,0]*self.pgas_vel[:,1]-self.pgas_coo[:,1]*self.pgas_vel[:,0])/(np.sqrt(self.pgas_coo[:,0]**2 + self.pgas_coo[:,1]**2)) #angular velocity using the velocity and position in cartesian coordinates

        # Keep cartesian coordinates (x and y part) for easy postage stamp plotting
        self.pgas_cart_coo = np.copy(self.pgas_coo[:,:2])
        self.pgas_cart_vel = np.copy(self.pgas_vel[:,:2])
        # Compute total velocity (amplitude)
        self.pgas_veltot   = np.sqrt(np.sum(self.pgas_vel**2, axis=1)) # (N,)
        
        self.pgas_vel[:,0] = v_r #change coordinates
        self.pgas_vel[:,1] = v_phi #change coordinates
        self.pgas_coo[:,0] = r #change coordinates
        self.pgas_coo[:,1] = phi #change coordinates


        #== Work on stellar particles ==
        self.pstar_coo=np.dot(A,self.pstar_coo.T).T  #change coordinates
        self.pstar_vel=np.dot(A,self.pstar_vel.T).T 
        
        #Changing coordinates to cylindrical polar       
        phi= np.arctan(self.pstar_coo[:,1]/self.pstar_coo[:,0]) #phi component 
        phi[(self.pstar_coo[:,0] < 0)] += np.pi # if x is negative, need to add pi
        r  = np.sqrt(self.pstar_coo[:,0]**2 + self.pstar_coo[:,1]**2) #r component from cartesian to cylindrical
        v_r=(self.pstar_coo[:,0]*self.pstar_vel[:,0]+self.pstar_coo[:,1]*self.pstar_vel[:,1])/(np.sqrt(self.pstar_coo[:,0]**2 + self.pstar_coo[:,1]**2)) #radial velocity using the velocity and position in cartesian coordinates 
        v_phi=(self.pstar_coo[:,0]*self.pstar_vel[:,1]-self.pstar_coo[:,1]*self.pstar_vel[:,0])/(np.sqrt(self.pstar_coo[:,0]**2 + self.pstar_coo[:,1]**2)) #angular velocity using the velocity and position in cartesian coordinates

        #_coo = np.copy(self.pgas_coo)#coordinates before the change in order to check everything functions well
        # Keep cartesian coordinates (x and y part) for easy postage stamp plotting
        self.pstar_cart_coo = np.copy(self.pstar_coo[:,:2])
        self.pstar_cart_vel = np.copy(self.pstar_vel[:,:2])
        # Compute total velocity (amplitude)
        self.pstar_veltot   = np.sqrt(np.sum(self.pstar_vel**2, axis=1)) # (N,)
        
        self.pstar_vel[:,0] = v_r #change coordinates
        self.pstar_vel[:,1] = v_phi #change coordinates
        self.pstar_coo[:,0] = r #change coordinates
        self.pstar_coo[:,1] = phi #change coordinates