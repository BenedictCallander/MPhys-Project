import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
import pwlf
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as pat

df =pd.read_csv("csv/alldatadone.csv")
snapshots = list(df['snapshot'])
inner = list(df['slope1'])
valids = [21,33,50,67,78,91,99]
maxvals = [] ; minvals = [] ; meanvals = []
dftemp = pd.read_csv("csv/alldatadone.csv")
dftemp = dftemp.groupby(['snapshot']).mean().reset_index()
for i in valids:
    dftry =pd.read_csv('csv/alldatadone2.csv')
    dftry= dftry[dftry['snapshot'].isin([i])]
    maxvals.append(max(dftry['slope2']))
    minvals.append(min(dftry['slope2']))
    meanvals.append(np.mean(dftry['slope2']))
print(maxvals)
print(minvals)
print(meanvals)


    
'''
plt.figure(figsize=(20,12))
plt.title("Evolution of Inner slope for TNG50-1 subhalos",fontsize=25)
plt.plot(valids,maxvals,'k+',ms=10)
plt.plot(dftemp['snapshot'],abs(dftemp['slope1']),'r--')
plt.xlabel("Snapshot",fontsize=20);plt.ylabel("Inner gradient",fontsize=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
plt.savefig('iner_Ev.png')
plt.close()
'''