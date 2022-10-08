import pandas as pd 
import requests


ids =[]; links = []; sfr =[]
url = []
mass_gas = []; mass = []; mass_stars = []
pos_x = []; pos_y = []; pos_z = []
starmetallicity = []
gasmetallicity = []
gasmetallicitysfrweighted = []
for i in range (5)
    result = r.get ()
    links.append(url)
    sfr.append(result['sfr'])
    mass_gas.append(result['mass_gas'])
    mass.append(result['mass'])
    mass_stars.append(result['mass_stars'])
    pos_x.append(result['pos_x'])
    pos_y.append(result['pos_y'])
    pos_z.append(result['pos_y'])
    starmetallicity.append(result['starmetallicity'])
    gasmetallicity.append(result['gasmetallicity'])
    gasmetallicitysfrweighted.append(result['gasmetallicitysfrweighted'])
df = pd.DataFrame({
    "id" : ids,
    "link" : url,
    "sfr" :  sfr,
    "mass_stars" : mass_stars,
    "mass_total" : mass,
    "mass_gas" :  mass_gas,
    "posx" :  pos_x,
    "posy" :  pos_y,
    "posz" :  pos_z,
    "metallicity_stars" : starmetallicity,
    "metallicity_gas" : gasmetallicity,
    "m_gas_sfr_weighted" : gasmetallicitysfrweighted})


