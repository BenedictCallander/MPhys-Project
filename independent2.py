#%%
#attempt to obtain data for multiple subhalos in z=0, tng 50
import numpy as np
import requests
import pandas as pd
import h5py
import matplotlib.pyplot as plt
baseurl = baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}

def get(path, params = 'None'):
    # function to make a HTTP GET request to a defined path
    r = requests.get(path, params = params, headers = headers)

    r.raise_for_status()
    

    if r.headers['content-type'] == 'application/json':
        return r.json()
    
    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string
    return r


base_url = 'https://www.tng-project.org/api/TNG50-1/snapshots/z=0/subhalos/'
sim_metadata = get(base_url)
cutout_request = {'gas':'gasmetallicitysfr'}

metallicity = []
masslog =[]
ilist=[]
for i in range(859077):
    fileurl = base_url+str(i)
    met = get(fileurl)
    metallicity.append(met['gasmetallicity'])
    masslog.append(met['mass_gas'])
    ilist.append(i)
    #print(i, met['starmetallicity'])
df = pd.DataFrame({'i': ilist, 'mass': masslog, 'metallicity': metallicity})
pd.DataFrame.to_csv('ind1.csv')

plt.figure()
plt.plot(masslog,metallicity, 'r.')
plt.xscale('log')
plt.yscale('log')
plt.show()
plt.close()


#starmetallicitymaxrad
#gasmetallicitymaxrad

# %%
