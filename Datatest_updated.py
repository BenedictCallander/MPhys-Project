#%% 
import numpy as np
import requests
import h5py
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

r = get(baseUrl)
r.keys()
['simulations']
len(r['simulations'])


names = [sim['name'] for sim in r['simulations']]

# print(names)


i = names.index('Illustris-3')

sim = get(r['simulations'][i]['url'])
sim.keys()

print(sim['num_dm'])

snaps = get(sim['snapshots'])
#print('length of snapshot', len(snaps)) #find total number of snapshots
#print(snaps[-1]) #print -1 index to find total number of snapshots

snap = get(snaps[-1]['url'])
print(snap)



subs = get(snap['subhalos'], {'limit':220})#obtain list of subhalos inside selected snapshot 
print(subs.keys(),'\n \n')
print(subs['results'][0])
print(subs['next'])
#test to see if github working
#black panther trailer


#request first 20 subhalos, sort by descending stellar mass



subs = get(snap['subhalos'],{'limit':20, 'order_by':'-mass_stars'})
#for i in range(5):
    #print(subs['results'][i]['id'])

sub = get(subs['results'][1]['url'])
#print(sub['cm_x'])

url = sub['related']['parent_halo']+'info.json'
#print(url)



#now request main progenitor branch from sublink merger trees of this subhalo

mpb1 = get(sub['trees']['sublink_mpb']) #file saved, mpb1 contrains filename

f = h5py.File(mpb1,'r')
print(f.keys())
print(len(f['SnapNum']))
print(f['SnapNum'][:])


mpb2 = get( sub['trees']['lhalotree_mpb'] ) # file saved, mpb2 contains the filename

with h5py.File(mpb2,'r') as f:
    print(len(f['SnapNum']))


with h5py.File(mpb2,'r') as f:
    pos = f['SubhaloPos'][:]
    snapnum= f['SnapNum'][:]
    subid = f['SubhaloNumber'][:]
#for i in range (3):
 #   plt.plot(snapnum,pos[:,i] - pos[0,i], label=['x','y','z'][i])
#plt.legend()
#plt.xlabel('Snapshot Number')
#plt.ylabel('Pos$_{x,y,z}$(z) - Pos(z=0)')

url = sim['snapshots']+'z=1'
print(url)

snap = get(url)

i = np.where(snapnum==85)
print(subid[i])

sub_prog_url = "http://www.tng-project.org/api/Illustris-3/snapshots/85/subhalos/185/"
sub_prog = get(sub_prog_url)

cutout_request = {'gas':'Coordinates,Masses'}
cutout = get(sub_prog_url+'cutout.hdf5', cutout_request)

with h5py.File(cutout,'r') as f:
    x = f['PartType0']['Coordinates'][:,0] - sub_prog['pos_x']
    y = f['PartType0']['Coordinates'][:,1] - sub_prog['pos_y']
    dens = np.log10(f['PartType0']['Masses'][:])

plt.hist2d(x,y,weights=dens,bins=[150,100])
plt.xlabel('$\Delta x$ [ckpc/h]')
plt.ylabel('$\Delta y$ [ckpc/h]');
# %%
