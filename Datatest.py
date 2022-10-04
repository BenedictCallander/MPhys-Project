#%% 

import requests

baseurl = baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}

def get(path, params = 'None'):
    # function to make a HTTP GET request to a defined path
    r = requests.get(path, params = params, headers = headers)

    r.raise_for_status()
    

    if r.headers['content-type'] == 'application/json':
        return r.json()
    return r

r = get(baseUrl)
r.keys()
['simulations']
len(r['simulations'])

r['simulations'][0]
{'name': 'Illustris-1',
 'num_snapshots': 134,
 'url': 'http://www.tng-project.org/api/Illustris-1/'}

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

# %%
