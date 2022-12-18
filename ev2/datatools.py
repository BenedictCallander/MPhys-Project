import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import scipy 



class evslopes:
    def __init__(self) -> None:
        pass
    def shapechar(x,y,z):
        '''
        take inputs of 3 slopes -> return value based on slope 
        '''
        
        #shape 1: steep negative inner followed by shallow negative mid/outer
        if x<0 and y<0 and abs(x)<abs(y):
            return 1
        #shape 2: shallow negative inner followed by steeper negative mid/outer
        elif x<0 and y<0 and abs(y)<abs(x):
            return 2
        #shape 4:positive inner followed by positive outer( Something has gone wrong)
        elif x>0 and y>0:
            return 4
        #shape 5: steep inner down followed by shallow negative outer:
        elif x<-0.5 and y<0 and y>-0.1:
            return 5
        #shape6: steep inner up followed by shallow negative outer:
        elif x>0.5 and y<0 and y>-0.1:
            return 6
        #shape 3: positive inner followed by negative outers:
        elif x>0 and y<0:
            return 3
        else:
            return 7
        
    def genchar(x,y):
        #shape 1: steep negative inner:
        if x<-0.3 and y<0:
            return 1 
        #shape 2: steep positive inner:
        elif x>0.3 and y<0:
            return 2
        #shape 3: anything else 
        else:
            return 3 
    def geninnerdir(x):
        if x>0:
            return 1
        elif x<0: 
            return 2
        else:
            return 3
df = pd.read_csv("csv/tng67MSslopes.csv")
df = df.dropna()
xs=list(df['slope1'])
ys=list(df['slope2'])
zs=list(df['slope3'])

values=[]
for i in range(len(xs)):
    value = evslopes.geninnerdir(xs[i])
    values.append(value)

df.insert(8,'Values',values)
df.sort_values(by='Values',inplace=True)
df.to_csv("try672.csv")

vals = [1,2,3]
for i in vals:
    df = pd.read_csv("try67.csv")
    df=df[df['Values'].isin([i])]
    d = list(df['Values'])
    print("number of subhalos with shape {} = {}".format(i,len(d)))
    
'''
    
'''

