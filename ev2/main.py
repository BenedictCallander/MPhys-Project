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

dfin = pd.read_csv("csv/alldatadonebreak.csv")

snapshots=[21,33,50,67,78,91,99]
'''
for i in snapshots:
    dfsnap = dfin[dfin['snapshot'].isin([i])].copy()
    x = list(dfsnap['breakpoint'])
    vals = []
    for j in range(len(x)):
        val = BCUTILS.shapeanalysis.get_inner_direction(x[j])
        vals.append(val)
    dfsnap.insert(7,"shape",vals,True)
    dfsnap.to_csv("shapes{}.csv".format(i)) 
    print("done for {}".format(i))

for i in snapshots:
    df = pd.read_csv("shapes{}.csv".format(i))
    df1 = df[df['shape'].isin([1])].copy()
    df2 = df[df['shape'].isin([2])].copy()
    df3 = df[df['shape'].isin([3])].copy()
    print("{} Total {}: number of:\n 1:{}\n 2:{}\n 3:{}".format(i,len(df.index),len(df1.index)/len(df.index),len(df2.index)/len(df.index),len(df3.index)/len(df.index)))

'''

for i in snapshots:
    dfsnap = dfin[dfin['snapshot'].isin([i])].copy()
    x = list(dfsnap['breakrad'])
    print("snapshot {} mean breakpoint is {}".format(i,np.mean(x)))
    #print("snapshot {} median breakpoint is {}".format(i,np.median(x)))

'''
sub = subhalo("TNG50-1",99,117250)
sub.galcen()
sub.ang_mom_align('gas')
sub.rad_transform()
dfg = sub.df_gen('gas','comb')
dfg2 = sub.combfilter(dfg,10)

print(dfg2)

subvis=plots(19)
subvis.met_histogram(dfg2,'Y')

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
'''