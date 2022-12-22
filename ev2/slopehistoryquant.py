import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pwlf

#hdf5 binary file manipulation
import h5py

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


def docombine(i):
    try:
        path1 = "files/historycutouts/evdir_{}/slope{}.csv".format(i,i)
        path2 = "files/historycutouts/evdir_{}/historydata_{}.csv".format(i,i)
        df1 = pd.read_csv(path1)
        df2 = pd.read_csv(path2)
        df2.insert(6,"slope1",df1['slope1'])
        df2.insert(7,"slope2",df1['slope2'])
        path = "files/historycutouts/evdir_{}/alldata{}.csv".format(i,i)
        df2.to_csv(path)
        print("done for {}".format(i))
    except FileNotFoundError as e:
        print("no file for {}".format(i))
        return print(e)

        
    
dfrun =pd.read_csv("csv/traceids.csv")
ids = list(dfrun['id'])

returns = Parallel(n_jobs=20)(delayed(docombine)(i) for i in ids)
dataframes = []
for i in ids:
    try:
        path="files/historycutouts/evdir_{}/alldata{}.csv".format(i,i)
        df = pd.read_csv(path)
        dataframes.append(df)
    except FileNotFoundError:
        print("1")
    
dfall = pd.concat(dataframes)
dfall.to_csv("alldatadone2.csv")