from itertools import groupby
import logging
from random import random
from re import sub # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il
from scipy.signal import medfilt
from scipy.optimize import curve_fit
from joblib import Parallel, delayed
import os
from scipy.signal import savgol_filter


df = pd.read_csv("test1.csv")

def line(x):
    y = pow(10,((2*x)-20.5))
    return y

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
        print("invalid")
print(len(valids))
df_in = pd.read_csv("test1.csv")
df2 = df[df['id'].isin(valids)]
df_in = df_in[df_in['sfr']>0]
valid_id = list(df_in['id'])
print(df2)
df2.to_csv("mainseq2.csv")
df3 = df_in[df_in['id'].isin(valids)]
df3.to_csv("testing2.csv")

#'''
xval = np.linspace(0,13,100)
def line(m,x,b):
    y = 10**((m*x)+b)
    return y 

df2=pd.read_csv("slopeplot.csv")
plt.figure(figsize=(20,12))
plt.plot(xval, line(2,xval,-20.5), 'r-', label = "y=$10^{mx+b}$")
#plt.scatter((df_in['mass']),(df_in['sfr']),c=df2['slope'],cmap='viridis',vmin=-0.2,label = 'Subhalos')
plt.scatter((df3['mass']), df3['sfr'],c = df2['slope'], cmap='viridis',label = 'subhalos_MS')
plt.xlabel("Subhalo Mass (log Msun)", fontsize=20)
plt.yscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15); plt.legend(loc='upper right')
plt.xlim(7.5,13)
plt.ylim(10e-6,10e2)
plt.ylabel("Log(Subhalo SFR)",fontsize=20)
plt.title("Galaxy Classification with Metallicity gradient visualisation",fontsize=20)
plt.colorbar().set_label(label = "Metallicity Linear Fit Slope",size=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.savefig("slope3.png")
plt.close()
#'''
'''
plt.figure(figsize=(20,12))
plt.plot(df3['mass'],df3['met'],'r+')
plt.savefig('met.png')
plt.close()
'''