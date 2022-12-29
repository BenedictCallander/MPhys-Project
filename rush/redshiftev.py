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

df33= pd.read_csv("csv/tng33bkpc.csv")
df67= pd.read_csv("csv/tng67bkpc.csv")
df99= pd.read_csv("csv/tng99bkpc.csv")

[7,8,9,10,11,12]
lower = 7
upper = 8
i=0

def do(dfin):
    df = dfin.copy()
    df = df.dropna()
    slope1 = list(df['slope1'])
    slope2 = list(df['slope2'])
    mean1 = np.var(slope1)
    mean2 = np.var(slope2)
    print("mean slope1: {} slope2: {}".format(mean1,mean2))

do(df33)
do(df67)
do(df99)