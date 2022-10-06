#%%
#
# Test code - filtering of only SFR >0 - reduce dataset tp 17553
#

from email.mime import base
import numpy as np 
import requests
import h5py 
import pandas as pd 
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

search_q = "?sfr__gt=0.0"

base_url = "http://www.tng-project.org/api/TNG50-2/snapshots/z=0/subhalos/"
url = base_url+search_q
valid_subs = get(url)
print(valid_subs['count'])
print(len(valid_subs['results']))
valid_ids = [ valid_subs['results'][i]['id'] for i in range(100)]
print(valid_ids)

mass = []
metallicity = []
for i in valid_ids:
    url = base_url+str(i)
    data = get(url)
    mass.append(data['mass_gas'])
    metallicity.append(data['gasmetallicity'])
    print("id # ", i, 'complete')

df = pd.DataFrame({'id': valid_ids, 'mass': mass, 'metallicity': metallicity})
df.to_csv('data1.csv')


    


# %%
