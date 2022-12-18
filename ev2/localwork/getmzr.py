#%%
import logging
from random import random # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
from joblib import Parallel, delayed
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il


headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
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


class onsub:
    def __init__(self,subID):
        self.subID = subID
        self.subURL = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/{}/".format(subID)
    
    def dataget(self):
        subhalo = get(self.subURL)
        mass = subhalo['mass_log_msun']
        met = subhalo['gasmetallicitysfrweighted']
        subID = subhalo['id']
        return (subID,mass,met)

dfin = pd.read_csv("tng99MAIN.csv")
ids = list(dfin['id'])
def dofunc(i):
    obj = onsub(i)
    subID,mass,met = obj.dataget()
    print("done for {}".format(i))
    return (subID,mass,met)

returns = Parallel(n_jobs=-1)(delayed(dofunc)(i) for i in ids)
df = pd.DataFrame(returns, columns = ['subID','mass','met'])
df.to_csv("MZR99.csv")