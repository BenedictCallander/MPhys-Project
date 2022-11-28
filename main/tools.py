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


class analyse:
    def __init__(self,subhalos, slopes, MainSeq):
        self.subDF = pd.read_csv(subhalos)
        self.slopeDF = pd.read_csv(slopes)
        self.mainDF = pd.read_csv(MainSeq)    
        
    def findX(df,lower,upper):
        r'''
        findX: Set X limits to query data in specific mass ranges 
        '''
        df = df.dropna()
        df = df[df['mass']>lower]
        df = df[df['mass']<upper]
        slopes = list((df['slope']))
        mean = np.median(slopes)
        return print("TNG33 MSQ in MASS Range {} -> {} median slope is :{}".format(lower,upper,mean))
    
    def findY(df,lower,upper):
        r'''
        findY: Set Y limits to query data in specific mass ranges 
        '''
        df = df.dropna()
        df = df[df['sfr']>lower]
        df = df[df['sfr']<upper]
        slopes = list((df['slope']))
        mean = np.median(slopes)
        return print("TNG33 MSQ in MASS Range {} -> {} median slope is :{}".format(lower,upper,mean))
    
    def boxfind(df, minX,maxX,minY,maxY):
        r'''
        boxfind: define box on M-SFR plot to select subalo subset and analyse slope statistics 
        
        INPUTS:
        
        df: pandas DataFrame
        
        contains subhalo data for plotting on M/SFR, including metallicity gradient data 
        
        minX: float
        
        position for lower X filter to be set 
        
        maxX: float 
        
        position for upper X filter to be set 
        
        minY: float
        
        position for lower Y filter to be set 
        
        maxY:
        
        position for upper Y filter to be set 

        
        
        '''
        df = df.dropna()
        df = df[df['mass'].between(minX,maxX)]
        df = df[df['sfr'].between(minY,maxY)]
        slopes = list(df['slope'])
        mean = np.mean(slopes)
        minslope = min(slopes)
        maxslope = max(slopes)
        sloperange = maxslope-minslope
        
        print("slope values for M {}->{}, SFR {}->{} \n".format(minX,maxX,minY,maxY))
        print("median slope in box = {} : minimum slope = {} : maximum slope = {}".format(mean, minslope,maxslope))
        print("range in slopes = {}".format(sloperange))
        return (slopes,mean,minslope,maxslope,sloperange)
    
    