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
print(sfrsubs['count'])
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


#df_analysis.to_csv('csv/tng33subhalos.csv')
xval = np.linspace(0,13,100)
yvals = np.linspace(10e-6,10e3,100)
def line(m,x,b):
    y = 10**((m*x)+b)
    return y 


plt.figure(figsize=(15,10))
plt.plot(mass_sfr,sfr_sfr,'g+')
#sns.kdeplot(x=mass_sfr, y=np.log10(sfr_sfr))
plt.plot(xval, line(2,xval,-20), 'r-', label = "y=$10^{mx+b}$")
plt.plot(xval, line(2,xval,-18.5), 'g-', label = "y=$10^{mx+b}$")
plt.axvline(x=8.5)
plt.axvline(x=9.5)
plt.axhline(y=10e-2)
plt.yscale('log')
plt.ylabel('log(SFR)')
plt.ylim(10e-6,10e2)
plt.xlim(7,13)
plt.xlabel('Mass (log10 Msun)')
plt.show()
plt.close()



def line1(x):
    y = pow(10,((2*x)-20.5))
    return y
def line2(x):
    y =pow(10,(2*x)-18.5)
    return y

df = df_analysis.copy()

df = df[df['mass']<11]
df = df[df['mass']>]
df = df[df['sfr']<1]
df = df[df['sfr']>0.2]
print(df)
df.to_csv('traceids.csv')

'''
def MSfilterup(dfin):
    df = dfin
    ids = list(df['id'])
    masses = list(df['mass'])
    sfr = list((df['sfr']))
    valids = []
    for i in range(len(ids)):
        value=line1((masses[i]))
        if value<((sfr[i])):
            valids.append(ids[i])
            print(i)
        else:
            continue
    return valids
valid1 = MSfilterup(df)
df = df[df['id'].isin(valid1)]


def MSfilterdown(dfin):
    df = dfin
    ids = list(df['id'])
    masses = list(df['mass'])
    sfr = list((df['sfr']))
    valids = []
    for i in range(len(ids)):
        value=line2((masses[i]))
        if value>((sfr[i])):
            valids.append(ids[i])
            print(i)
        else:
            continue
    return valids

valid2 = MSfilterdown(df)
df = df[df['id'].isin(valid2)]
df.to_csv('ids99.csv')
plt.figure(figsize=(15,10))
plt.plot(mass_sfr,sfr_sfr,'g+')
#sns.kdeplot(x=mass_sfr, y=np.log10(sfr_sfr))
plt.plot(xval, line1(xval), 'r-', label = "y=$10^{mx+b}$")
plt.plot(xval, line2(xval), 'g-', label = "y=$10^{mx+b}$")
plt.plot(df['mass'],df['sfr'],'r*')
plt.axvline(x=8.5)
plt.axvline(x=9.5)
plt.yscale('log')
plt.ylabel('log(SFR)')
plt.ylim(10e-6,10e2)
plt.xlim(7,13)
plt.savefig("classif.png")
plt.xlabel('Mass (log10 Msun)')
plt.close()
df.to_csv("treefinds.csv")
print(list(df['id']))
end = time.time()
print("runtime :{}".format(end-start))
'''