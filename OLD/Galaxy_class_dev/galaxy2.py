#%%
#
# Test code - filtering of only SFR >0 - reduce dataset tp 17553
#
import logging # http logging for debugging purpouses
import time #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 
import csv

class galaxy:
    def __init__(self, num,id,link,sfr,mass_stars,mass_total,mass_gas,posx,posy,posz,metallicity_stars,metallicity_gas,m_gas_sfr_weighted):
        self.num=num
        self.id=id
        self.link=str(link)
        self.sfr=sfr
        self.mass_stars=mass_stars
        self.mass_total=mass_total
        self.mass_gas =mass_gas
        self.posx= posx
        self.posy=posy
        self.posz=posz
        self.metallicity_stars=metallicity_stars
        self.metallicity_gas=metallicity_gas
        self.m_gas_sfr_weighted= m_gas_sfr_weighted
    def pos_3D(self):
        pos3d = [self.posx,self.posy,self.posz] 
        return pos3d

galaxy1 = galaxy(9,37,'http://www.tng-project.org/api/TNG50-2/snapshots/z=0/subhalos/37',0.757937,0.139988,6.74045,0.37798,6917.02,24525.6,24525.6,0.011673,0.00921499,0.016004)
print(galaxy1.pos_3D())


with open('file.csv', newline = '') as csv_file:
    reader = csv.reader(csv_file)
    next(reader,None)
    for 