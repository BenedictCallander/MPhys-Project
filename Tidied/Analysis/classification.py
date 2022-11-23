#%%
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


baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/33/subhalos/"
sfr_q = "?limit=77655&sfr__gt=0.001"


sfrurl = baseurl+sfr_q
sfrsubs = get(sfrurl)
mass=[]
sfr=[]

sfr_ids = [sfrsubs['results'][i]['id'] for i in range(sfrsubs['count'])]
mass_sfr = []
sfr_sfr = []
urls = []
for i,id in enumerate(sfr_ids):
    mass_sfr.append(sfrsubs['results'][i]['mass_log_msun'])
    sfr_sfr.append(sfrsubs['results'][i]['sfr'])
    urls.append(sfrsubs['results'][i]['url'])

df_analysis = pd.DataFrame({
    "id": sfr_ids,
    "mass": mass_sfr,
    "sfr": sfr_sfr,
    "url": urls
})
df_analysis.to_csv('csv/tng33subhalos.csv')
'''
xval = np.linspace(0,13,100)
def line(m,x,b):
    y = 10**((m*x)+b)
    return y 
def line2(m,x,b):
    y = (m*x) + b
    return y 

plt.figure(figsize=(15,10))
plt.plot(mass_sfr,np.log10(sfr_sfr),'g+')
#sns.kdeplot(x=mass_sfr, y=np.log10(sfr_sfr))
plt.plot(xval, line2(2,xval,-20), 'r-', label = "y=$10^{mx+b}$")
plt.yscale('log')
plt.ylabel('log(SFR)')
plt.ylim(-6.5,2)
plt.xlim(7,14)
plt.xlabel('Mass (log10 Msun)')
plt.savefig('png/classification/SFR_M_TNG50-1_991.png')
plt.close()
'''
end = time.time()
print("runtime :{}".format(end-start))
# %%
