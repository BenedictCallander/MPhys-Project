import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pwlf
import scipy.stats as stats
#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il

#Own module containing utility functions 

#runtime calculation 
import time

#Computational functions - simultaneous calculations to make use of multi-core CPU
from joblib import Parallel, delayed

#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter

df33= pd.read_csv("csv/ffs33.csv")
df67= pd.read_csv("csv/ffs67.csv")
df99= pd.read_csv("csv/ffs99.csv")


def do(dfin):
    lower = 7
    upper = 8
    slopes = [] ; vars = []
    while upper<=12:
        df = dfin.copy()
        df = df.dropna()
        df.mass = np.log10(df.mass/0.7)
        df = df[df['mass']<upper].copy()
        df = df[df['mass']>lower].copy()
        slope = list(df['slope'])
        mean = np.mean(slope)
        slopes.append(mean)
        vars.append(np.var(slope)) 
        print("from mass {} slope = {}: var = {}: Nslopes = {}".format(lower, mean,np.var(slope),len(slope)))
        lower = upper
        upper = upper+1
do(df99)