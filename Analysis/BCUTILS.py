import time  # runtime calculation import numpy as np #data handling
from random import random

import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests  # obtain data from API server
import seaborn as sns 
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()

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

def line(x):
    y = pow(10,((2*x)-20.5))
    return y

def MSfilter(dfin,dfcombin, saveloc):
    '''
    INPUTS:
        *dfin -> main dataframe containing all subhalos, mass and SFR (API DATA)
        *dfcombin -> secondary dataframe from galaxy processing script, with slope column
        *saveloc -> filepath for output file to be saved to
        
    PROCESSES: 
        *Filtering -> uses MS line equation and performs calculations to determine whether subhalo lies above or below line, appends IDS of 
                      valid subhalos to a list
        *Data Restructuring -> uses list of Main Sequence subhalo IDS to filter dataframes such that only valid subhalos remain, inserts columns for slope
                        and AIC value for these dataframes and saves output to file 
    OUTPUTS:
        * Data for Main Sequence subhalos saved in .csv format to location specified by saveloc
        * Upon completion saved file path is printed 
    '''
    df = dfin
    ids = list(df['id'])
    masses = list(df['mass'])
    sfr = list((df['sfr']))
    valids = []
    valids_m = []
    valid_sfr = []
    for i in range(len(ids)):
        value=line((masses[i]))
        if value<((sfr[i])):
            valids.append(ids[i])
            print(i)
        else:
            continue
    df_in =df
    dfslope = dfcombin
    dfslope = dfslope[dfslope['id'].isin(valids)]
    df3 = df_in[df_in['id'].isin(valids)]
    df3.insert(2,"slope",dfslope['slope'],True)
    df3.insert(3,"AIC", dfslope['AICval'],True)
    df3.insert(4,"met", dfslope['met'],True)
    df3.to_csv(saveloc)
    return print("dataframe saved to: ", saveloc)

def subhalo_classification(snapshot,contours):
    if snapshot == 33:
        baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/33/subhalos/"
        search_q = "?limit=77655&sfr__gt=0.0"
        url = baseurl+ search_q
    elif snapshot == 67:
        baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/67/subhalos/"
        search_q = "?limit=26980&sfr__gt=0.0"
        url = baseurl+ search_q
    elif snapshot ==99:
        baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/"
        search_q = "?limit=17553&sfr__gt=0.0"
        url = baseurl+ search_q
    xvals = np.linspace(0,13,1000)
    
    subs = get(url)
    ids = [subs['results'][i]['id'] for i in range(subs['count'])]
    mass = []
    sfr = []
    urls = []
    for i,id in enumerate(ids):
        mass.append(subs['results'][i]['mass_log_msun'])
        sfr.append(subs['results'][i]['sfr'])
        urls.append(subs['results'][i]['url'])
    
    plt.figure(figsize=(20,12))
    plt.plot(mass,np.log10(sfr),'g+')
    if contours=='Y':
        sns.kdeplot(x=mass, y=np.log10(sfr))
    elif contours =='N':
        print('No Contours selected')
    plt.plot(xvals, line(2,xvals,-20.5), 'r-', label = "y=$10^{mx+b}$")
    #plt.yscale('log')
    plt.ylabel('log(SFR)')
    plt.ylim(-6.5,2)
    plt.xlim(7,14)
    plt.xlabel('Mass (log10 Msun)')
    plt.savefig('png/classification/SFR_M_TNG50-1_{}.png'.format(snapshot))
    plt.close()
dfin = pd.read_csv("csv/tng99subhalos.csv")
df2=pd.read_csv("csv/tng99slopes.csv")
MSfilter(dfin,df2,'csv/tng99MAIN.csv')