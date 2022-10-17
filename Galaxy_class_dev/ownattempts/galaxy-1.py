import numpy as np 
import requests
import h5py 
import pandas as pd 
import scipy 
import matplotlib.pyplot as plt


def get(path, params=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

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