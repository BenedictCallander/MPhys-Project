#%%
#
# Test code - filtering of only SFR >0 - reduce dataset tp 17553
#
import logging # http logging for debugging purpouses
import time #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 

baseurl = baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
#logging.basicConfig(level=logging.DEBUG)
r = requests.Session() #- think ineffective? 
start = time.time() # start runtime calculator for code improvement 

#define utility function - to get information from TNG API 
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

#define search criteria for subhalo filtering 
search_q = "?limit=200&sfr__gt=0.0"

#set base url as subhalos list 
base_url = "http://www.tng-project.org/api/TNG50-2/snapshots/z=0/subhalos/"
url = base_url+search_q
#find all subhalos that fill search criteria 
valid_subs = get(url)


#print number of valid subhalos 
print("number of star forming subhalos in simulation = ",valid_subs['count'])
#print(len(valid_subs['results']))

#find the subhalo ID numbers of each of these star forming subhalos 
valid_ids = [ valid_subs['results'][i]['id'] for i in range(100)]

#
#-currently limited to 1 page, 100 results per page, could test changing limits for runspeed
#

#
# read url for 'next' page from API, load next page and collect ID numbers of valid subhalos
#

for i in range(20): #range = number of pages past 1 to include in data download
    urlnext = valid_subs['next']
    valid_subs= get(urlnext)
    ids2 = [ valid_subs['results'][i]['id'] for i in range(100)]
    valid_ids.extend(ids2) #add ID numbers to collection
    ids2.clear()

ids =[]; links = []; sfr =[]
url = []
mass_gas = []; mass = []; mass_stars = []
pos_x = []; pos_y = []; pos_z = []
starmetallicity = []
gasmetallicity = []
gasmetallicitysfrweighted = []
ticker = 1

print("sending API request for ", len(valid_ids), "subhalos")
for i in valid_ids:
    url=base_url+str(i)
    #print(url)
    result = get(url)
    links.append(url)
    sfr.append(result['sfr'])
    mass_gas.append(result['mass_gas'])
    mass.append(result['mass'])
    mass_stars.append(result['mass_stars'])
    pos_x.append(result['pos_x'])
    pos_y.append(result['pos_y'])
    pos_z.append(result['pos_y'])
    starmetallicity.append(result['starmetallicity'])
    gasmetallicity.append(result['gasmetallicity'])
    gasmetallicitysfrweighted.append(result['gasmetallicitysfrweighted'])
    ticker = ticker + 1
    pc = (ticker/(len(valid_ids)))*100
    if pc%5 ==0:
        print('Request {:.2f} % complete, {:.2f} s'.format(pc, (time.time()-start)))
    if pc ==5:
        print('projected runtime = {:.2f}'.format(pc*20))

df = pd.DataFrame({
    "id" : valid_ids,
    "link" : links,
    "sfr" :  sfr,
    "mass_stars" : mass_stars,
    "mass_total" : mass,
    "mass_gas" :  mass_gas,
    "posx" :  pos_x,
    "posy" :  pos_y,
    "posz" :  pos_z,
    "metallicity_stars" : starmetallicity,
    "metallicity_gas" : gasmetallicity,
    "m_gas_sfr_weighted" : gasmetallicitysfrweighted})
df.to_csv('bigdata1.csv')

end = time.time()
print('runtime is', (end-start),' seconds')
# %%
