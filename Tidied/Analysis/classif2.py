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

df33 = pd.read_csv('csv/tng33subhalos.csv')
df67 = pd.read_csv('csv/tng67subhalos.csv')
df99 = pd.read_csv('csv/tng99subhalos.csv')

mass_33 = list(df33['mass'])
mass_67 = list(df67['mass'])
mass_99 = list(df99['mass'])

sfr_33 = list(df33['sfr'])
sfr_67 = list(df67['sfr'])
sfr_99 = list(df99['sfr'])
xval = np.linspace(0,13,100)

def line(m,x,b):
    y = 10**((m*x)+b)
    return y 

def line2(m,x,b):
    y = (m*x) + b
    return y 

def plotgraph(mass, sfr,i):
    plt.figure(figsize=(15,10))
    plt.plot(mass,np.log10(sfr),'g+')
    sns.kdeplot(x=mass, y=np.log10(sfr),   )
    plt.plot(xval, line2(2,xval,-20), 'r-', label = "y=$10^{mx+b}$")
    #plt.yscale('log')
    plt.ylabel('log(SFR)')
    plt.ylim(-6.5,2)
    plt.xlim(7,14)
    plt.xlabel('Mass (log10 Msun)')
    filename = 'png/classification/SFR_M_TNG50-1_{}.png'.format(i)
    plt.savefig(filename)
    plt.close()