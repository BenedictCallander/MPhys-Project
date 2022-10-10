#%%
#
# Test code - filtering of only SFR >0 - reduce dataset tp 17553
#
import logging #http debugging
import time #runtime calculations
import numpy as np #data management 
import requests #API services
import h5py  #large scale binary file format management
import pandas as pd #data management
import matplotlib.pyplot as plt #data plotting

baseurl = baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}

logging.basicConfig(level=logging.DEBUG) #debugging configuration
r = requests.Session() #open continuous connection to baseURL
r.get(baseurl)

start = time.time()
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
print("number of star forming subhalos in simulation = ",valid_subs['count'])
#print(len(valid_subs['results']))
valid_ids = [ valid_subs['results'][i]['id'] for i in range(100)]
for i in range(2):
    urlnext = valid_subs['next']
    valid_subs= get(urlnext)
    ids2 = [ valid_subs['results'][i]['id'] for i in range(100)]
    valid_ids.extend(ids2)#
    ids2.clear()
#print(urlnext)
print(valid_ids)
end = time.time()
print("Request runtime - {:.4f}s".format(end-start))
base_url = "http://www.tng-project.org/api/TNG50-2"
params = {'gas':'Masses,GFM_Metallicity'}
for i in valid_ids:
    file_url = base_url +"files/snapshot-99" + str(i)+'.hdf5'
    saved_filename = get(file_url, params)
    print(saved_filename)

    

# %%
