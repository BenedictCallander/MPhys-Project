import os 
import time  # runtime calculation import numpy as np #data handling
from random import random
from re import sub  # http logging for debugging purpouses
import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests  # obtain data from API server
import h5py
from joblib import Parallel,delayed

headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()
pd.options.mode.chained_assignment = None  # default='warn'

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



baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/?order_by=-mass_stars"
subs = get(baseurl)
subIDS = [subs['results'][i]['id'] for i in range(12,25)]

for i,id in enumerate(subIDS):
    url = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/{}/".format(i)
    sub = get(url)
    mpb1 = get(sub['trees']['sublink_mpb'] )
    mpb2 = get( sub['trees']['lhalotree_mpb'] )
    print("files downloaded for {}".format(i))
    print("runtime {}".format((time.time())-start))

