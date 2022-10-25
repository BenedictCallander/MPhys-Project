import logging
from random import random # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
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

baseurl = "https://www.tng-project.org/api/TNG100-1/snapshots/70/subhalos/"
sq = '?sfr__gt=0.0'
sfr_q = "?limit=105181&sfr__gt=0.0"
all_q = "?limit=4716081"


sfrurl = baseurl+sfr_q
sfrsubs = get(sfrurl)
allurl = baseurl+all_q
allsubs = get(allurl)
'''
allsub = get(baseurl)
print(allsub['count']) 
>>4716081

print(sfrsubs['count'])
>>105181
all_ids = [ allsubs['results'][i]['id'] for i in range(4716081)]

for i,id in enumerate(all_ids):
    mass.append(allsubs['results'][i]['mass_log_msun'])
    sfr.append(allsubs['results'][i]['sfr'])
'''
mass=[]
sfr=[]

sfr_ids = [ sfrsubs['results'][i]['id'] for i in range(105181)]
mass_sfr = []
sfr_sfr = []
for i,id in enumerate(sfr_ids):
    mass_sfr.append(sfrsubs['results'][i]['mass_log_msun'])
    sfr_sfr.append(sfrsubs['results'][i]['sfr'])
    
plt.figure(figsize=(15,10))
plt.plot(sfr_sfr, mass_sfr,'b.')
plt.xlabel('SFR')
plt.ylabel('Mass (log10 Msun)')
plt.savefig('SFR_M_TNG100-1.png')
plt.close()

end = time.time()
print("runtime :{}".format(end-start))