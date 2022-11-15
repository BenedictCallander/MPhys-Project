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


df = pd.read_csv("csv/tng33subhalos.csv")

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
df_in =df
dfslope = pd.read_csv('csv/tng33slopes.csv')
dfslope = dfslope[dfslope['id'].isin(valids)]
df3 = df_in[df_in['id'].isin(valids)]
df3.insert(2,"slope",dfslope['slope'],True)
df3.to_csv("csv/tng33MSQ.csv")




'''
plt.figure(figsize=(20,12))
plt.plot(df3['mass'],df3['met'],'r+')
plt.savefig('met.png')
plt.close()
'''