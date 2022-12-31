import numpy as np 

import pandas as pd 


df33 = pd.read_csv("csv/tng33slopes.csv")
df67 = pd.read_csv("csv/tng67slopes.csv")
df99 = pd.read_csv("csv/tng99slopes.csv")

def do(dfin):
    df = dfin.copy()
    df=df.dropna()
    slope1 = list(df['slope'])
    print(np.mean(slope1))
    print(np.var(slope1))
    '''
    slope2 = list(df['slope2'])
    print(np.mean(slope2))
    print(np.var(slope2))
    '''

do(df33)

'''
lower = 7
upper = 8
slopes = [] ; vars = []
while upper<=12:
    df = pd.read_csv("csv/tng99slopes.csv")
    df = df.dropna()
    df.mass = np.log10(10e10*df.mass/0.7)
    df = df[df['mass']<upper].copy()
    df = df[df['mass']>lower].copy()
    slope = list(df['slope'])
    mean = np.mean(slope)
    slopes.append(mean)
    vars.append(np.var(slope)) 
    print("from mass {} slope = {}: var = {}: Nslopes = {}".format(lower, mean,np.var(slope),len(slope)))
    lower = upper
    upper = upper+1
'''

#print(slopes)
#print(vars)