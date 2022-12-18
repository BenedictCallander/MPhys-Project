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
xval = np.linspace(0,13,100)
yvals = np.linspace(10e-6,10e3,100)
def line(m,x,b):
    y = 10**((m*x)+b)
    return y 
#df_analysis.to_csv('csv/tng33subhalos.csv')
def line2(m,x,b):
    y = (m*x)+b
    return y

plt.figure(figsize=(15,10))
plt.plot(mass_sfr,np.log10(sfr_sfr),'g+',label = "TNG50-1 snapshot 99 Subhalos")
sns.kdeplot(x=mass_sfr, y=np.log10(sfr_sfr))
plt.plot(xval, line2(2,xval,-20.5), 'r-', label = "y=2$M_\odot$-20.5")
#plt.plot(xval, line(2,xval,-18.5), 'g-', label = "y=$10^{mx+b}$")
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.xlabel("Total Mass [$M_\odot$]",fontsize=20)
plt.ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)

plt.ylim(-6,3)
plt.xlim(7,14)
plt.legend(loc='upper right')
plt.savefig("classif2.png")
plt.close()

end=time.time()
print("runtime is {}".format(end-start))