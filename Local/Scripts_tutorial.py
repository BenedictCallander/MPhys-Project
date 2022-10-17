import logging # http logging for debugging purpouses
import time #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 
import illustris_python as il

saved_filename =r'C:\Users\Administrator\Desktop\MPhys Project\Local\offsets_099.hdf5'
with h5py.File(saved_filename) as f:
    # NOTE! If the subhalo is near the edge of the box, you must take the periodic boundary into account! (we ignore it here)
    dx = f['PartType0']['Coordinates'][:,0] - sub['pos_x']
    dy = f['PartType0']['Coordinates'][:,1] - sub['pos_y']
    dz = f['PartType0']['Coordinates'][:,2] - sub['pos_z']
    electronabnd = f['PartType0']['ElectronAbundance'][:]
    nhydrogenabnd = f['PartType0']['NeutralHydrogenAbundance'][:]
    density = f['PartType0']['Density'][:]
    gfmmetallicity = f['PartType0']['GFM_Metallicity'][:]
    intenergy = f['PartType0']['InternalEnergy'][:]
    sfr = f['PartType0']['StarFormationRate'][:]
    mass = f['PartType0']['Masses'][:]