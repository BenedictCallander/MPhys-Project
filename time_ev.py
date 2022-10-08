#%%

import random
import logging # http logging for debugging purpouses
import time #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
#import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 

start = time.time()


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
base_url = "http://www.tng-project.org/api/TNG50-2/snapshots/z=0/subhalos/"

sq = '?sfr_gt=0.0'
url1 = base_url+sq
valid = get(url1)
length = valid['count']
print(length)
search_q = "?limit=200&sfr__gt=0.0"
url = base_url+search_q
valid_subs = get(url)
print(url)
valid_ids = [ valid_subs['results'][i]['id'] for i in range(10057)]
print(len(valid_ids))

mass = []
sfr = []

for i,id in enumerate(valid_ids):
    mass.append(valid_subs['results'][i]['mass_log_msun'])
    sfr.append(valid_subs['results'][i]['sfr'])

df = pd.DataFrame({"id": valid_ids, "sfr": sfr, "mass": mass})

df.to_csv('test1.csv')
colors = [plt.cm.hsv(i) for i in np.linspace(0, 1, 600)]

plt.figure(figsize = (20,12), dpi = 500)
plt.scatter(df['mass'], df['sfr'], c = df['id'], cmap = 'hsv')
plt.yscale('log')
plt.ylabel("log(SFR)", fontsize = 25)
plt.xlabel("Mass (Log Msun)", fontsize = 25)
plt.tick_params(axis = 'both', which = 'both',direction = 'inout', length = 15, labelsize = 15)
plt.colorbar(label='arbitrary 3rd')
plt.savefig('sfr.png')
plt.show()
plt.close()

# %%
