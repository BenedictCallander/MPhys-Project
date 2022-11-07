#%%

import time #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
#import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 
from joblib import Parallel,delayed
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
#'''
base_url = "http://www.tng-project.org/api/TNG100-1/snapshots/70/subhalos/"
massive_ids = []

search_q = "?limit=105181&sfr__gt=0.0"
vurl = base_url+search_q
valid_subs = get(vurl)
print(vurl)
valid_ids = [ valid_subs['results'][i]['id'] for i in range(105181)]
valid_urls = []

idvals=[]
urlvals =[]
for i,id in enumerate(valid_ids):
    idvals.append(valid_subs['results'][i]['id'])
    urlvals.append(valid_subs['results'][i]['url'])
    
df=pd.DataFrame({"id": idvals, "url":urlvals})
df.to_csv("parallel.csv")
#'''



start = time.time()
df = pd.read_csv("parallel.csv")
def gfm_get

def gfm_get (url):
    try:
        subdata = get(url)
        ids=(subdata['id'])
        mass=(subdata['mass_gas'])
        GFM_met=(subdata['gasmetallicity'])
        radius = subdata['halfmassrad']
        sfr = subdata['sfr']
        print(ids)
    except:
        print("Connection refused by the server..")
        print("Let me sleep for 5 seconds")
        print("ZZzzzz...")
        time.sleep(5)
        print("Was a nice sleep, now let me continue...")
    return ids, GFM_met, mass,radius,sfr

returns = Parallel(n_jobs=20)(delayed(gfm_get)(url) for url in df["url"])
print(len(returns))
print(np.ndim(returns))
print(returns)
df2=pd.DataFrame(returns,columns=['ids','met','mass','radius','sfr'])
print(df2)
df2.to_csv('rad.csv')
end=time.time()
print("parallel time = {}".format(end-start))


'''
start2=time.time
ids = []
GFM_met = []
mass = []
radius = []
ticker = 0
total = len(df['url'])
print(total)
for i in df['url']:
    subdata= get(i)
    ids.append(subdata['id'])
    mass.append(subdata['mass_gas'])
    GFM_met.append(subdata['gasmetallicity'])
    radius.append(subdata['halfmassrad'])
    print("Subhalo number {} \n".format((ticker/total))*100)
    ticker=ticker+1

df3 = pd.DataFrame({"id":ids, "met":GFM_met, "mass":mass})
df3.to_csv("rad.csv")


end2 = time.time()
print ("runtime = {}".format(end2-start2))
'''
# %%
