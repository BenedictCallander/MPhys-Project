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

df = pd.read_csv('mainseq2.csv')
df2 = pd.read_csv('testing2.csv')

plt.figure(figsize=(20,12))
plt.scatter((df2['mass']),(12+np.log10(df['met'])), c = (df['slope']), cmap = 'viridis', vmin=-0.15, label = 'Main Sequence Subhalos')
plt.xlabel("Subhalo Mass (log Msun)", fontsize=20)
plt.ylabel("12+ $log_{10}$ ${O}/{H}$",fontsize=20)
#plt.xscale('log')
plt.ylim(7.5,12)
plt.title("Mass Metallicity Relation for Main Sequence Subhalos:",fontsize=20)
plt.colorbar().set_label(label = "Metallicity Linear Fit Slope",size=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15); plt.legend(loc='upper right')
plt.savefig('met.png')
plt.close()