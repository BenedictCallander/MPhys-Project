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
df = pd.read_csv('csv/tng33subhalos.csv')
df1 = pd.read_csv('csv/tng33slopes.csv')

xvals = np.linspace(0,13,100)

def line(m,x,b):
    y = pow(10,((m*x)+b))
    return y 

print(np.mean(df1['slope']))
plt.figure(figsize=(20,12))
plt.scatter((df['mass']),df['sfr'], c = (df1['slope']), cmap = 'viridis',vmin=-0.25, label = 'Main Sequence Subhalos')
plt.plot(xvals,line(2,xvals,-20.5),'r-')
plt.xlabel("Subhalo Mass (log Msun)", fontsize=20)
plt.ylabel("log(SFR)")
#plt.ylabel("12+ $log_{10}$ ${O}/{H}$",fontsize=20)
plt.yscale('log')
plt.ylim(10e-6, 10e2)
plt.xlim(7.5,12)
plt.title("Mass Metallicity Relation for Main Sequence Subhalos:",fontsize=20)
plt.colorbar().set_label(label = "Metallicity Linear Fit Slope",size=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15); plt.legend(loc='upper right')
filename = 'png/slopegrads/slope_67.png'
plt.savefig(filename)
plt.close()