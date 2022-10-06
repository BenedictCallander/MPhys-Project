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
search_q = "?sfr__gt=0.0"

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

for i in range(1):
    urlnext = valid_subs['next']
    valid_subs= get(urlnext)
    ids2 = [ valid_subs['results'][i]['id'] for i in range(100)]
    valid_ids.extend(ids2) #add ID numbers to collection
    ids2.clear()
    
#debugging print statements 
#print(urlnext)
#print(valid_ids)

#
#UI message to broadcast program state 
#
print("sending API Fetch request for", len(valid_ids),"datapoints")
mass = []
metallicity = []
ticker = 1
for i in valid_ids:
    url = base_url+str(i)
    data = get(url)
    mass.append(data['mass_gas'])
    metallicity.append(data['gasmetallicity'])
    ticker = ticker+1
    pc = (ticker/len(valid_ids))*100
    #
    #Progress update broadcast
    #UI message to broadcast program state 
    #only display 5% increments to stop visual overload
    #
    if pc%5 ==0:
        print("Request is {:.2f} % complete".format(pc))

#
#concentrate data into dataframe for printing to  and then 
#
df = pd.DataFrame({'id': valid_ids, 'mass': mass, 'metallicity': metallicity})
df.to_csv('data_large.csv')
end = time.time()
print("request complete - data written to .csv file")
print("Request runtime - {:.4f}s".format(end-start))

plt.figure()
plt.plot(mass,metallicity, 'g.')
plt.xscale('log')
plt.yscale('log')
plt.show()
plt.close()



# %%
