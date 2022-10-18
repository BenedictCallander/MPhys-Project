from random import random
import numpy as np 
import matplotlib.pyplot as plt 
import requests
import logging 
import scipy 
import pandas as pd 
import illustris_python as il 
import h5py

#
#BCUTILS
#

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

def subhalo_filter():
    #search_q = "?sfr__gt=0.0"

#set base url as subhalos list 
    base_url = "http://www.tng-project.org/api/TNG50-2/snapshots/z=0/subhalos/?limit=1000&sfr__gt=0.0"
    url = base_url#+search_q
    #find all subhalos that fill search criteria 
    valid_subs = get(url)


    #print number of valid subhalos 
    print("number of star forming subhalos in simulation = ",valid_subs['count'])
    #print(len(valid_subs['results']))

    #find the subhalo ID numbers of each of these star forming subhalos 
    valid_ids = [ valid_subs['results'][i]['id'] for i in range(100)]
    return valid_ids
valid_ids=subhalo_filter()
print(valid_ids)

def visualise_cutout(id, type, lim):


    sub_prog_url = "http://www.tng-project.org/api/TNG50-2/snapshots/99/subhalos/"+str(id)+"/"
    sub_prog = get(sub_prog_url)
    
    cutout_request = {'gas':'Coordinates,Masses,GFM_Metallicity'}
    cutout = get(sub_prog_url+"cutout.hdf5", cutout_request)


    with h5py.File(cutout,'r') as f:
        x = f['PartType0']['Coordinates'][:,0] - sub_prog['pos_x']
        y = f['PartType0']['Coordinates'][:,1] - sub_prog['pos_y']
        if (type=='gas'):
            dens =(f['PartType0']['Masses'][:])
        elif(type =='met'): 
            dens =(f['PartType0']['GFM_Metallicity'][:])

    plt.figure()
    plt.hist2d(x,y,weights=dens,bins=[1000,500], cmap = 'inferno', vmin = min(dens), vmax = max(dens))
    plt.xlabel('$\Delta x$ [ckpc/h]')
    plt.ylabel('$\Delta y$ [ckpc/h]')
    plt.xlim(-lim,lim)
    plt.ylim(-lim,lim)
    plt.savefig('hist_met_{}.png'.format(id))
    plt.close()

    return print("graph plotted")
random_ids = [9938,9681,9513,7401,5676,5511,6317]

visualise_cutout(8,'gas',10)

