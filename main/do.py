import numpy as np
import pandas as pd 



df33= pd.read_csv("csv/tng33MSslopes.csv")
df67= pd.read_csv("csv/tng67MSslopes.csv")
df99= pd.read_csv("csv/tng99MSslopes.csv")

[7,8,9,10,11,12]
lower = 7
upper = 8
i=0

'''
while upper<=12:
    df = df99.copy()
    df = df.dropna()
    df = df[df['mass']>lower]
    df = df[df['mass']<upper]
    slopes = list(df['slope'])
    print("mean for mass min {} is {}: N slopes is {}".format(lower, np.mean(slopes),len(slopes)))
    lower = lower+1
    upper=upper+1
'''


def do(dfin):
    df = dfin.copy()
    df = df.dropna()
    slope1 = list(df['inside'])
    slope2 = list(df['outside'])
    mean1 = np.mean(slope1)
    mean2 = np.mean(slope2)
    print("mean slope1: {} slope2: {}".format(mean1,mean2))

do(df33)
do(df67)
do(df99)