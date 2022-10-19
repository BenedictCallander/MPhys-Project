#
# BCUTILS Function Library for MPhys Project 
#
import logging # http logging for debugging purpouses
import time #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 
import random
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}


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


def subhalo_filter(limit):
    #search_q = "?sfr__gt=0.0"

#set base url as subhalos list 
    base_url = "http://www.tng-project.org/api/TNG50-2/snapshots/z=0/subhalos/?limit={}&sfr__gt=0.0".format(limit)
    url = base_url#+search_q
    #find all subhalos that fill search criteria 
    valid_subs = get(url)


    #print number of valid subhalos 
    print("number of star forming subhalos in simulation = ",valid_subs['count'])
    #print(len(valid_subs['results']))

    #find the subhalo ID numbers of each of these star forming subhalos 
    valid_ids = [ valid_subs['results'][i]['id'] for i in range(limit)]
    return valid_ids


def visualise_cutout(id, type, lim):


    sub_prog_url = "http://www.tng-project.org/api/TNG50-2/snapshots/99/subhalos/"+str(id)+"/"
    sub_prog = get(sub_prog_url)
    
    cutout_request = {'gas':'Coordinates,Masses,GFM_Metallicity'}
    cutout = get(sub_prog_url+"cutout.hdf5", cutout_request)


    with h5py.File(cutout,'r') as f:
        x = f['PartType0']['Coordinates'][:,0] - sub_prog['pos_x']
        y = f['PartType0']['Coordinates'][:,1] - sub_prog['pos_y']
        if (type=='gas'):
            dens =np.exp(f['PartType0']['Masses'][:])**4
        elif(type =='met'): 
            dens =-np.log10(f['PartType0']['GFM_Metallicity'][:])**3

    plt.figure()
    plt.scatter(x,y,c=-np.log10(dens),cmap='inferno',s=1.2)
    #plt.hist2d(x,y,weights=dens,bins=[1500,1000], cmap = 'afmhot', vmin = min(dens), vmax = (max(dens)))
    plt.xlabel('$\Delta x$ [ckpc/h]')
    plt.ylabel('$\Delta y$ [ckpc/h]')
    #plt.xlim(-10,10)
    #plt.ylim(-7,4)
    plt.savefig('hist_met_{}_{}.png'.format(type,id))
    plt.close()

    return print("graph plotted for subhalo{}".format(id))

visualise_cutout(40,'met',5)