#
# This wasnt a good idea 
#

import requests
import h5py 
import numpy as np 
import matplotlib.pyplot as plt 
import scipy
import pandas as pd 
import time 
from joblib import Parallel, delayed
from scipy.signal import medfilt
from scipy.optimize import curve_fit 

headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
baseurl='https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/'
def get(path, params = None):
    '''
    API request function - expansion upon requests.get() to provide error codes and .hdf5 file reading 
    Inputs:
        - Path - URL for request to fetch data from 
        - Params - additional parameters that can be used (i.e. search and limit in the case of the ILLUSTRIS TNG API  )
    '''
    #Make API request for path (needs to be str variable) including header (API KEY) and other params
    r = requests.get(path, params=params,headers=headers) 

    #check for HTTP error codes
    r.raise_for_status()

    #detect content type for function to read 

    if r.headers['content-type'] == 'application/json':
        return r.json()
    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

class APIutils:
    def __init__(self, url,):
        self.url=url
    def masssfrplot(self,criteria, filename):
        baseurl = self.url
        url = baseurl+str(criteria)
        subsdata=get(url)
        ids = []
        mass=[]
        sfr=[]
        valid_ids = [ subsdata['results'][i]['id'] for i in range(17553)]
        for i,id in enumerate(valid_ids):
            mass.append(subsdata['results'][i]['mass_log_msun'])
            ids.append(subsdata['results'][i]['id'])
            sfr.append(subsdata['results'][i]['sfr'])

        plt.figure(figsize=(15,10))
        plt.plot(mass, np.log10(sfr), 'k+')
        plt.xlabel("Mass (Log Msun)")
        plt.ylabel("log(SFR)")
        plt.savefig(filename)
        plt.close()


plot1=APIutils(baseurl)    
plot1.masssfrplot("?limit=17553&sfr__gt=0.0", "testing.png")
