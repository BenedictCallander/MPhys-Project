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
        
    def find(df,lower,upper):
        df = df.dropna()
        df = df[df['mass']>lower]
        df = df[df['mass']<upper]
        slopes = list((df['slope']))
        mean = np.median(slopes)
        return print("TNG33 MSQ in MASS Range {} -> {} median slope is :{}".format(lower,upper,mean))