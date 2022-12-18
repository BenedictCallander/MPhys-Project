import logging
from random import random # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il
import seaborn as sns 

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


baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/"
sfr_q = "?limit=17553&sfr__gt=0.0"
sfrurl = baseurl+sfr_q
sfrsubs = get(sfrurl)
mass=[]
sfr=[]
print(sfrsubs['count'])
sfr_ids = [sfrsubs['results'][i]['id'] for i in range(sfrsubs['count'])]
mass_sfr = []
sfr_sfr = []
urls = []
for i,id in enumerate(sfr_ids):
    mass_sfr.append(sfrsubs['results'][i]['mass_log_msun'])
    sfr_sfr.append(sfrsubs['results'][i]['sfr'])
    urls.append(sfrsubs['results'][i]['url'])

df = pd.DataFrame({
    "id": sfr_ids,
    "mass": mass_sfr,
    "sfr": sfr_sfr,
    "url": urls
})
df.to_csv('tng99subhalos.csv')

plt.figure(figsize=(20,12))
plt.plot( df['mass'],(df['sfr']), 'g+')
plt.yscale('log')
plt.xlabel("Total Mass [$M_\odot$]",fontsize=20);plt.ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(9.5,11.5)
plt.savefig("classif.png")