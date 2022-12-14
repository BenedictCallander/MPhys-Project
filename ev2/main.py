import matplotlib.pyplot as plt; import numpy as np ;import pandas as pd
import pwlf ;import h5py
#Read data from web API and monitor HTTP traffic 
import requests  
#specialised functions to query Illustris TNG data 
import illustris_python as il
#specific functions for fitting utilities
from scipy.optimize import curve_fit ; from scipy.signal import medfilt, savgol_filter

import BCUTILS 
from subobj import subhalo

from BCUTILS import UTILITY,plots
'''
sub = subhalo("TNG50-1",99,19)
sub.galcen()
sub.ang_mom_align('gas')
sub.rad_transform()
dfg = sub.df_gen('gas','comb')
dfg2 = sub.combfilter(dfg,10)

print(dfg2)

subvis=plots(19)
subvis.gas_visual(dfg2,4)
'''

dfin = pd.read_csv("csv/tng67MAIN.csv")
ids = list(dfin['id'])
print(len(ids))
dataframes=[]
invalids = []
valids = []
for i in ids:
    try:
        fpath = "dfdump/df{}.csv".format(i)
        df=pd.read_csv(fpath)
        dataframes.append(df)
        valids.append(i)
    except FileNotFoundError:
        invalids.append(i)
print(len(invalids))    
data2 = pd.concat(dataframes)
data2.to_csv("all67working.csv")

dfin=dfin[dfin['id'].isin(valids)]
dfin.to_csv("csv/tng67MAIN.csv")