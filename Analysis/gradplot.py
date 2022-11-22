
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


MS1 = pd.read_csv('csv/tng33MSQ.csv')
MS2 = pd.read_csv('csv/tng67MSQ.csv')
MS3 = pd.read_csv('csv/tng99MSQ.csv')
df1.slope = (10*((df1.slope-df1.slope.min())/(df1.slope.max()-df1.slope.min())))
df3.slope = (10*((df3.slope-df3.slope.min())/(df3.slope.max()-df3.slope.min())))
df5.slope = (10*((df5.slope-df5.slope.min())/(df5.slope.max()-df5.slope.min())))

MS1.slope = 10*((MS1.slope-MS1.slope.min())/(MS1.slope.max()-MS1.slope.min()))
MS2.slope = 10*((MS2.slope-MS2.slope.min())/(MS2.slope.max()-MS2.slope.min()))
MS3.slope = 10*((MS3.slope-MS3.slope.min())/(MS3.slope.max()-MS3.slope.min()))

xvals = np.linspace(0,13,100)

def line(m,x,b):
    y = pow(10,((m*x)+b))
    return y 
'''
plt.figure(figsize=(20,12))
plt.scatter((df['mass']),df['sfr'], c = (df1['slope']), cmap = 'magma',vmin=5,vmax=9, label = 'Main Sequence Subhalos')
plt.plot(xvals,line(2,xvals,-20.5),'r-')
plt.xlabel("Subhalo Mass (log Msun)", fontsize=20)
plt.ylabel("log(SFR)", fontsize =20)
#plt.ylabel("12+ $log_{10}$ ${O}/{H}$",fontsize=20)
plt.yscale('log')
plt.ylim(10e-4, 10e2)
plt.xlim(7.5,12)
plt.title("Subhalo Classification snapshot 99",fontsize=20)
plt.colorbar().set_label(label = "Metallicity Linear Fit Slope",size=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15); plt.legend(loc='upper right')
filename = 'png/slopegrads/met99.png'
plt.savefig(filename)
plt.close()

'''


font = {'family': 'serif',
        'weight': 'normal',
        'size': 16,
        }

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
fig,((ax0,ax1,ax2),(ax3,ax4,ax5),(ax6,ax7,ax8)) = plt.subplots(nrows = 3, ncols = 3,figsize=(30,24))
divider0 = make_axes_locatable(ax0); divider1 = make_axes_locatable(ax1);divider2 = make_axes_locatable(ax2);divider3 = make_axes_locatable(ax3)
divider4 = make_axes_locatable(ax4);divider5 = make_axes_locatable(ax5);divider6 = make_axes_locatable(ax6)
divider7 = make_axes_locatable(ax7);divider8 = make_axes_locatable(ax8)
cax0 = divider0.append_axes('right', size='5%', pad=0.05)
cax1 = divider1.append_axes('right', size='5%', pad=0.05)
cax2 = divider2.append_axes('right', size='5%', pad=0.05)
cax3 = divider3.append_axes('right', size='5%', pad=0.05)
cax4 = divider4.append_axes('right', size='5%', pad=0.05)
cax5 = divider5.append_axes('right', size='5%', pad=0.05)
cax6 = divider6.append_axes('right', size='5%', pad=0.05)
cax7 = divider7.append_axes('right', size='5%', pad=0.05)
cax8 = divider8.append_axes('right', size='5%', pad=0.05)

ax0.set_yscale('log')
ax0.set_title("Snap_033: Z=2")
im0 = ax0.scatter((df['mass']),df['sfr'], c = (df1['slope']), cmap = 'magma',label = 'Main Sequence Subhalos')
ax0.plot(xvals,line(2,xvals,-20.5),'r-')
ax0.set_ylim(10e-6, 10e2)
ax0.set_xlim(7.5,12)
ax0.set_ylabel("log(SFR) $log_{10}$", math_fontfamily='stix')


ax1.set_yscale('log')
ax1.set_title("Snap_067: Z=0.5")
im1= ax1.scatter((df2['mass']),df2['sfr'], c = (df3['slope']), cmap = 'magma',label = 'Main Sequence Subhalos')
ax1.plot(xvals,line(2,xvals,-20.5),'r-')
ax1.set_ylim(10e-6, 10e2)
ax1.set_xlim(7.5,12)

ax2.set_yscale('log')
ax2.set_title("Snap_099: Z=0")
im2 = ax2.scatter((df4['mass']),df4['sfr'], c = (df5['slope']), cmap = 'magma', label = 'Main Sequence Subhalos')
ax2.plot(xvals,line(2,xvals,-20.5),'r-')
ax2.set_ylim(10e-6, 10e2)
ax2.set_xlim(7.5,12)

ax3.set_ylabel("log(GFM_Metallicity)")
#ax3.set_yscale('log')
ax3.set_title("Metallicity 33: Z=2")
im3 = ax3.scatter((df['mass']),12+np.log10(dfalt1['met']), c = (df1['slope']), cmap = 'magma', label = 'Main Sequence Subhalos')
#ax3.set_ylim(10e-8, 10e2)
ax3.set_xlim(7.5,12)

#ax4.set_yscale('log')
ax4.set_title("Metallicity 67: Z=0.5")
im4 = ax4.scatter((df2['mass']),12+np.log10(dfalt3['met']), c = (df3['slope']), cmap = 'magma', label = 'Main Sequence Subhalos')
#ax4.set_ylim(10e-8, 10e2)
ax4.set_xlim(7.5,12)

#ax5.set_yscale('log')
ax5.set_title(" Metallicity 99: Z=0 ")
im5 = ax5.scatter((df4['mass']),12+np.log10(dfalt5['met']), c = (df5['slope']), cmap = 'magma', label = 'Main Sequence Subhalos')
#ax5.set_ylim(10e-8, 10e2)
ax5.set_xlim(7.5,12)

ax6.set_yscale('log')
ax6.set_title("Snap_033: Z=2")
im6 = ax6.scatter((MS1['mass']),MS1['sfr'], c = MS1['slope'], cmap = 'magma',vmin=2, label = 'Main Sequence Subhalos')
ax6.plot(xvals,line(2,xvals,-20.5),'r-')
ax6.set_ylim(10e-6, 10e2)
ax6.set_xlim(7.5,12)
ax6.set_ylabel("log(SFR)")

ax7.set_yscale('log')
ax7.set_title("Snap_033: Z=2")
im7 = ax7.scatter((MS2['mass']),MS2['sfr'], c = MS2['slope'], cmap = 'magma', label = 'Main Sequence Subhalos')
ax7.plot(xvals,line(2,xvals,-20.5),'r-')
ax7.set_ylim(10e-6, 10e2)
ax7.set_xlim(7.5,12)
ax7.set_ylabel("log(SFR)")

ax8.set_yscale('log')
ax8.set_title("Snap_033: Z=2")
im8 = ax8.scatter((MS3['mass']),MS3['sfr'], c =MS3['slope'], cmap = 'magma', label = 'Main Sequence Subhalos')
ax8.plot(xvals,line(2,xvals,-20.5),'r-')
ax8.set_ylim(10e-6, 10e2)
ax8.set_xlim(7.5,12)
ax8.set_ylabel("log(SFR)")


fig.colorbar(im0, cax=cax0, orientation='vertical')
fig.colorbar(im1, cax=cax1, orientation='vertical')
fig.colorbar(im2, cax=cax2, orientation='vertical')
fig.colorbar(im3, cax=cax3, orientation='vertical')
fig.colorbar(im4, cax=cax4, orientation='vertical')
fig.colorbar(im5, cax=cax5, orientation='vertical')
fig.colorbar(im6, cax=cax6, orientation='vertical')
fig.colorbar(im7, cax=cax7, orientation='vertical')
fig.colorbar(im8, cax=cax8, orientation='vertical')

fig.tight_layout()
fig.subplots_adjust(top=0.9)
fig.suptitle("Redshift progression of galaxy classification and metallicity slope", fontsize=20)
fig.savefig('png/slopegrads/combined2.png')

