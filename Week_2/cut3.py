import logging # http logging for debugging purpouses
import time #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 
baseurl = baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start = time.time()
def get(path, params=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r
id=0

sub_prog_url = "http://www.tng-project.org/api/TNG100-2/snapshots/99/subhalos/"+str(id)+"/"
sub_prog = get(sub_prog_url)

cutout_request = {'gas':'Coordinates,Masses,GFM_Metallicity'}
cutout = get(sub_prog_url+"cutout.hdf5", cutout_request)


with h5py.File(cutout,'r') as f:
    x = f['PartType0']['Coordinates'][:,0] - sub_prog['pos_x']
    y = f['PartType0']['Coordinates'][:,1] - sub_prog['pos_y']
    dens =-np.log10(f['PartType0']['GFM_Metallicity'][:])
lim = 150
plt.figure()
plt.hist2d(x,y,weights=dens,bins=[5000,5000], cmap = 'inferno')
plt.xlabel('$\Delta x$ [ckpc/h]')
plt.ylabel('$\Delta y$ [ckpc/h]')
plt.xlim(-lim,lim)
plt.ylim(-lim,lim)
plt.savefig('hist_met_{}.png'.format(id))
plt.close()

end = time.time()
print('runtime {}'.format(end-start))
