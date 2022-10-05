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


base_url = "http://www.tng-project.org/api/Illustris-1/"
sim_metadata = get(base_url)
params = {'stars':'Coordinates,Masses,GFM_Metallicity'}

for i in range(sim_metadata['num_files_snapshot']):
    file_url = base_url + "files/snapshot-135." + str(i) + ".hdf5"
    saved_filename = get(file_url, params)
    print(saved_filename)
# %%
