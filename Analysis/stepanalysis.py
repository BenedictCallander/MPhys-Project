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

MS1df = pd.read_csv('csv/tng33MAIN.csv')
MS2df = pd.read_csv('csv/tng67MAIN.csv')
MS3df = pd.read_csv('csv/tng99MAIN.csv')
print(min(MS1df['slope']))
print(max(MS1df['slope']))
MS1df.slope = 10*((MS1df.slope-MS1df.slope.min())/(MS1df.slope.max()-MS1df.slope.min()))
MS2df.slope = 10*((MS2df.slope-MS2df.slope.min())/(MS2df.slope.max()-MS2df.slope.min()))
MS3df.slope = 10*((MS3df.slope-MS3df.slope.min())/(MS3df.slope.max()-MS3df.slope.min()))
MS1df['met'] = 12+np.log10(MS1df['met'] )
MS2df['met'] = 12+np.log10(MS2df['met'] )
MS3df['met'] = 12+np.log10(MS3df['met'] )

low = 7.5
up = 9.5
def find(df,lower,upper):
    df = df.dropna()
    df = df[df['mass']>lower]
    df = df[df['mass']<upper]
    slopes = list((df['slope']))
    mean = np.median(slopes)
    return print("TNG33 MSQ in MASS Range {} -> {} median slope is :{}".format(lower,upper,mean))
i=0
while i<3:
    find(MS1df,low,up)
    #find(MS2df,low,up)
    #find(MS3df,low,up)
    low = up;up=up+1
    i=i+1

    