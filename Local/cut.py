#
# Load progenitor subhalos to subset of TNG99 subhalos and compute their metallicity parameters
#

import numpy as np
import pandas as pd

#hdf5 binary file manipulation
import h5py

#Read data from web API and monitor HTTP traffic 
import requests  



#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}


fpath = "cutout_8.hdf5"
with h5py.File(fpath,'r') as f:
    print(f.keys())
