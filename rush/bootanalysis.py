import numpy as np 
import pandas as pd 

df = pd.read_csv("csv2/boot99.csv")


minvals = list(df['slope1'])
maxvals = list(df['slope2'])
slopes = list(df['slope'])

ranges = []
for i in range(len(minvals)):
    range = maxvals[i] - minvals[i]
    ranges.append(range)
    
print(np.mean(maxvals))
#print(np.mean(ranges)/np.mean(slopes))