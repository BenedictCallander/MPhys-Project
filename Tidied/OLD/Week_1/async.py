import asyncio
import os 
import time
import requests
import aiohttp

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

search_q = "?sfr__gt=0.0"

base_url = "http://www.tng-project.org/api/TNG50-2/snapshots/z=0/subhalos/"
url = base_url+search_q
valid_subs = get(url)
print("number of star forming subhalos in simulation = ",valid_subs['count'])
#print(len(valid_subs['results']))
valid_ids = [ valid_subs['results'][i]['id'] for i in range(100)]
for i in range(1):
    urlnext = valid_subs['next']
    valid_subs= get(urlnext)
    ids2 = [ valid_subs['results'][i]['id'] for i in range(100)]
    valid_ids.extend(ids2)#
    ids2.clear()
results = []

async def get_data():
    async with aiohttp.ClientSession as session:
        for i in valid_ids:
            #get request
            #get url
            #get subhalo data
            #write data to loops
            data = await session.get(url, ssl=False)
            results.append(await data.json())
    return data 

asyncio.run(get_data())

