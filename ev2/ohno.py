import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import pwlf
import seaborn as sns 
#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il

#Own module containing utility functions 
import BCUTILS

#runtime calculation 
import time

#Computational functions - simultaneous calculations to make use of multi-core CPU
from joblib import Parallel, delayed

#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter

Msun = 1.98847e30
#set basic constants during initialisation for easy 
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}


def get(path, params=None):
# make HTTP GET request to path
    headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
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

baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/?limit=17553&sfr__gt=0.0"

sfrsubs = get(baseurl)
idvals =[sfrsubs['results'][i]['id'] for i in range(sfrsubs['count'])]

def get_data(i):
    url = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/{}".format(i)
    sub = get(url)
    mass_half = sub['SubhaloMassInHalfRad']
    sfr_half = sub['sfrinhalfrad']
    met_half = sub['gasmetallicityhalfrad']
    subID = sub['id']
    print("got for subhalo {}".format(i))
    return (subID,met_half,sfr_half,mass_half)

returns = Parallel(n_jobs=30)(delayed(get_data)(i) for i in idvals)
df = pd.DataFrame(returns,columns = ['id','met','sfr','mass'])
df.to_csv("tng99.csv")