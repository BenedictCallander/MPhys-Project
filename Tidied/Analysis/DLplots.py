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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import savgol_filter
import matplotlib as mpl

df = pd.read_csv('csv/tng33subhalos.csv')
df1 = pd.read_csv('csv/tng33slopes.csv')
dfalt1= pd.read_csv('csv/tng33slopes.csv')
df2 = pd.read_csv('csv/tng67subhalos.csv')
df3 = pd.read_csv('csv/tng67slopes.csv')
dfalt3 = pd.read_csv('csv/tng67slopes.csv')

df4 = pd.read_csv('csv/tng99subhalos.csv')
df5 = pd.read_csv('csv/tng99slopes.csv')
df5 = pd.read_csv('csv/tng99slopes.csv')
dfalt5 = pd.read_csv('csv/tng99slopes.csv')


MS1 = pd.read_csv('csv/tng33MAIN.csv')
MS2 = pd.read_csv('csv/tng67MAIN.csv')
MS3 = pd.read_csv('csv/tng99MAIN.csv')
df1.slope = 10*((df1.slope-df1.slope.min())/(df1.slope.max()-df1.slope.min()))
df3.slope = 10*((df3.slope-df3.slope.min())/(df3.slope.max()-df3.slope.min()))
df5.slope = 10*((df5.slope-df5.slope.min())/(df5.slope.max()-df5.slope.min()))

MS1.slope = 10*((MS1.slope-MS1.slope.min())/(MS1.slope.max()-MS1.slope.min()))
MS2.slope = 10*((MS2.slope-MS2.slope.min())/(MS2.slope.max()-MS2.slope.min()))
MS3.slope = 10*((MS3.slope-MS3.slope.min())/(MS3.slope.max()-MS3.slope.min()))

Num1 = len(MS1['id'])
Num2 = len(MS2['id'])
Num3 = len(MS3['id'])

xvals = np.linspace(0,13,100)

def line(m,x,b):
    y = pow(10,((m*x)+b))
    return y 


plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

fig,(ax3,ax4,ax5) = plt.subplots(nrows = 1, ncols = 3,figsize=(30,8))
divider0 = make_axes_locatable(ax3); divider1 = make_axes_locatable(ax4);divider2 = make_axes_locatable(ax5)
cax0 = divider0.append_axes('right', size='5%', pad=0.05)
cax1 = divider1.append_axes('right', size='5%', pad=0.05)
cax2 = divider2.append_axes('right', size='5%', pad=0.05)


ax3.set_ylabel("12+$log_{10}(O/H)$")
ax3.set_xlabel("Mass $log_{10} M_{sun}$")
#ax3.set_yscale('log')
ax3.set_title("Metallicity 33: Z=2")
im3 = ax3.scatter((df['mass']),12+np.log10(dfalt1['met']), c = (df1['slope']), cmap = 'magma',vmin=5,vmax=7, label = 'Main Sequence Subhalos')
ax3.set_ylim(7,11.5)
ax3.text(10,8,"Number of Subhalos = {}".format(Num1),fontsize=15)
ax3.set_xlim(8,11.5)

#ax4.set_yscale('log')
ax4.set_title("Metallicity 67: Z=0.5")
im4 = ax4.scatter((df2['mass']),12+np.log10(dfalt3['met']), c = (df3['slope']), cmap = 'magma',vmin=6, label = 'Main Sequence Subhalos')
ax4.set_ylabel("12+$log_{10}(O/H)$")
ax4.set_xlabel("Mass $log_{10} M_{sun}$")
ax4.set_ylim(7,12)
ax4.text(10,8,"Number of Subhalos = {}".format(Num2),fontsize=15)
ax4.set_xlim(8,11.5)

#ax5.set_yscale('log')
ax5.set_ylabel("12+$log_{10}(O/H)$")
ax5.set_xlabel("Mass $log_{10} M_{sun}$")
ax5.set_title(" Metallicity 99: Z=0 ")
im5 = ax5.scatter((df4['mass']),12+np.log10(dfalt5['met']), c = (df5['slope']), cmap = 'magma',vmin=6, label = 'Main Sequence Subhalos')
ax5.set_ylim(7,12)
ax5.text(10,8,"Number of Subhalos = {}".format(Num3),fontsize=15)
ax5.set_xlim(8,11.5)



fig.colorbar(im3, cax=cax0, orientation='vertical')
fig.colorbar(im4, cax=cax1, orientation='vertical')
fig.colorbar(im5, cax=cax2, orientation='vertical')

fig.tight_layout()
fig.subplots_adjust(top=0.9, right = 0.9)
fig.suptitle("Redshift progression of TNG50 Subhalo Mass-Metallicity relation", fontsize=20)
fig.savefig('png/slopegrads/MASS_met_2_COMB.png')