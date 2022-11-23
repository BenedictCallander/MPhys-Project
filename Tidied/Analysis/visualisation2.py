import logging
import os
import time  # runtime calculation import numpy as np #data handling
from itertools import groupby
from random import random
from re import sub  # http logging for debugging purpouses

import h5py  # binary file manipulation
import illustris_python as il
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests  # obtain data from API server
from joblib import Parallel, delayed
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter

def line(m,x,b):
    y = (m*x)+b
    return y 

xvals = np.linspace(0,13,100)
    
def nineplot(saveloc):
    df = pd.read_csv('csv/tng33subhalos.csv')
    df1 = pd.read_csv('csv/tng33slopes.csv')

    df2 = pd.read_csv('csv/tng67subhalos.csv')
    df3 = pd.read_csv('csv/tng67slopes.csv')

    df4 = pd.read_csv('csv/tng99subhalos.csv')
    df5 = pd.read_csv('csv/tng99slopes.csv')

    MS1 = pd.read_csv('csv/tng33MSQ.csv')
    MS2 = pd.read_csv('csv/tng67MSQ.csv')
    MS3 = pd.read_csv('csv/tng99MSQ.csv')
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
    im0 = ax0.scatter((df['mass']),df['sfr'], c = (df1['slope']), cmap = 'viridis', vmin=-0.25,vmax = 0.05, label = 'Main Sequence Subhalos')
    ax0.plot(xvals,line(2,xvals,-20.5),'r-')
    ax0.set_ylim(10e-6, 10e2)
    ax0.set_xlim(7.5,12)
    ax0.set_ylabel("log(SFR)")

    ax1.set_yscale('log')
    ax1.set_title("Snap_067: Z=0.5")
    im1= ax1.scatter((df2['mass']),df2['sfr'], c = (df3['slope']), cmap = 'viridis', vmin=-0.2,label = 'Main Sequence Subhalos')
    ax1.plot(xvals,line(2,xvals,-20.5),'r-')
    ax1.set_ylim(10e-6, 10e2)
    ax1.set_xlim(7.5,12)

    ax2.set_yscale('log')
    ax2.set_title("Snap_099: Z=0")
    im2 = ax2.scatter((df4['mass']),df4['sfr'], c = (df5['slope']), cmap = 'viridis', vmin=-0.2, label = 'Main Sequence Subhalos')
    ax2.plot(xvals,line(2,xvals,-20.5),'r-')
    ax2.set_ylim(10e-6, 10e2)
    ax2.set_xlim(7.5,12)

    ax3.set_ylabel("log(GFM_Metallicity)")
    ax3.set_yscale('log')
    ax3.set_title("Metallicity 33: Z=2")
    im3 = ax3.scatter((df['mass']),df1['met'], c = (df1['slope']), cmap = 'viridis',vmax = 0.15, label = 'Main Sequence Subhalos')
    ax3.set_ylim(10e-8, 10e2)
    ax3.set_xlim(7.5,12)

    ax4.set_yscale('log')
    ax4.set_title("Metallicity 67: Z=0.5")
    im4 = ax4.scatter((df2['mass']),df3['met'], c = (df3['slope']), cmap = 'viridis', vmin=-0.2, label = 'Main Sequence Subhalos')
    ax4.set_ylim(10e-8, 10e2)
    ax4.set_xlim(7.5,12)

    ax5.set_yscale('log')
    ax5.set_title(" Metallicity 99: Z=0 ")
    im5 = ax5.scatter((df4['mass']),df5['met'], c = (df5['slope']), cmap = 'viridis', vmin=-0.2, label = 'Main Sequence Subhalos')
    ax5.set_ylim(10e-8, 10e2)
    ax5.set_xlim(7.5,12)

    ax6.set_yscale('log')
    ax6.set_title("Snap_033: Z=2")
    im6 = ax6.scatter((MS1['mass']),MS1['sfr'], c = (MS1['slope']), cmap = 'viridis',vmax = 0.1, label = 'Main Sequence Subhalos')
    ax6.plot(xvals,line(2,xvals,-20.5),'r-')
    ax6.set_ylim(10e-6, 10e2)
    ax6.set_xlim(7.5,12)
    ax6.set_ylabel("log(SFR)")

    ax7.set_yscale('log')
    ax7.set_title("Snap_033: Z=2")
    im7 = ax7.scatter((MS2['mass']),MS2['sfr'], c = (MS2['slope']), cmap = 'viridis',  label = 'Main Sequence Subhalos')
    ax7.plot(xvals,line(2,xvals,-20.5),'r-')
    ax7.set_ylim(10e-6, 10e2)
    ax7.set_xlim(7.5,12)
    ax7.set_ylabel("log(SFR)")

    ax8.set_yscale('log')
    ax8.set_title("Snap_033: Z=2")
    im8 = ax8.scatter((MS3['mass']),MS3['sfr'], c = (MS3['slope']), cmap = 'viridis', label = 'Main Sequence Subhalos')
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
    #fig.savefig('png/slopegrads/combined.png')
    fig.savefig(saveloc)
    return print(saveloc)
   
def tripleplot(saveloc):
    #dataframe inputs -> look for more efficient method 
    df = pd.read_csv('csv/tng33subhalos.csv')
    df1 = pd.read_csv('csv/tng33slopes.csv')

    df2 = pd.read_csv('csv/tng67subhalos.csv')
    df3 = pd.read_csv('csv/tng67slopes.csv')

    df4 = pd.read_csv('csv/tng99subhalos.csv')
    df5 = pd.read_csv('csv/tng99slopes.csv')
    #initialise suplot object and create 3 axes (ax0,ax1,ax2)
    fig,(ax0,ax1,ax2) = plt.subplots(nrows=1,ncols=3, figsize=(30,8))
    #make axes object locatable for overplotting of colorbars
    divider0 = make_axes_locatable(ax0); divider1 = make_axes_locatable(ax1);divider2 = make_axes_locatable(ax2)
    cax0 = divider0.append_axes('right', size='5%', pad=0.05)
    cax1 = divider1.append_axes('right', size='5%', pad=0.05)
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)

    ax0.set_yscale('log')
    ax0.set_title("Snap_033: Z=2")
    im0 = ax0.scatter((df['mass']),df['sfr'], c = (df1['slope']), cmap = 'viridis', vmin=-0.25,vmax = 0.05, label = 'Main Sequence Subhalos')
    ax0.plot(xvals,line(2,xvals,-20.5),'r-')
    ax0.set_ylim(10e-6, 10e2)
    ax0.set_xlim(7.5,12)
    ax0.set_ylabel("log(SFR)")

    ax1.set_yscale('log')
    ax1.set_title("Snap_067: Z=0.5")
    im1= ax1.scatter((df2['mass']),df2['sfr'], c = (df3['slope']), cmap = 'viridis', vmin=-0.2,label = 'Main Sequence Subhalos')
    ax1.plot(xvals,line(2,xvals,-20.5),'r-')
    ax1.set_ylim(10e-6, 10e2)
    ax1.set_xlim(7.5,12)

    ax2.set_yscale('log')
    ax2.set_title("Snap_099: Z=0")
    im2 = ax2.scatter((df4['mass']),df4['sfr'], c = (df5['slope']), cmap = 'viridis', vmin=-0.2, label = 'Main Sequence Subhalos')
    ax2.plot(xvals,line(2,xvals,-20.5),'r-')
    ax2.set_ylim(10e-6, 10e2)
    ax2.set_xlim(7.5,12)
    
    #plot customisation
    fig.colorbar(im0, cax=cax0, orientation='vertical')
    fig.colorbar(im1, cax=cax1, orientation='vertical')
    fig.colorbar(im2, cax=cax2, orientation='vertical')
    fig.tight_layout()
    fig.subplots_adjust(top=0.9)
    fig.suptitle("Redshift progression of galaxy classification and metallicity slope", fontsize=20)
    #fig.savefig('png/slopegrads/combined.png')
    fig.savefig(saveloc)
    return print(saveloc)