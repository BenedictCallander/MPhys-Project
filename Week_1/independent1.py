#%%
#attempt to obtain data for multiple subhalos in z=0, tng 50
import numpy as np
import requests
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


base_url = 'https://www.tng-project.org/api/TNG50-2/99/subhalos/'
sim_metadata = get(base_url)


metallicity = []
for i in range(859077):
    fileurl = baseurl+str(i)+'/cutout.hdf5'
    met = get(fileurl)
    met = met.json()
    print(met['gasmetallicitysfr'])

    


#starmetallicitymaxrad
#gasmetallicitymaxrad

# %%
