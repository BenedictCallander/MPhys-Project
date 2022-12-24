import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
import pwlf
from mpl_toolkits.axes_grid1 import make_axes_locatable
#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il

#runtime calculation 
import time


#Computational functions - simultaneous calculations to make use of multi-core CPU
from joblib import Parallel, delayed

#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()


'''
df1 = pd.read_csv("csv/tng33KPCslopesboth.csv")
df2= pd.read_csv("csv/tng67KPCslopesboth.csv")
df3 = pd.read_csv("csv/tng99KPCslopesboth.csv")


df3 = pd.read_csv("tng99medfilt.csv")

'''

df3 = pd.read_csv("tng99AIC2.csv")

df1 = df3[df3['AIC'].isin([1])].copy() #linear better
df2 = df3[df3['AIC'].isin([2])].copy() # piecewise better

plt.figure(figsize=(20,12),dpi=500)
plt.plot(df1['mass'],df1['sfr'],'k+', label = "Linear Fit",zorder=1)
plt.plot(df2['mass'],df2['sfr'],'g+', label = "Piecewise Fit",zorder=0)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xlabel("Log(Total Mass) [$M_\odot$]",fontsize=20)
plt.ylabel("Log(Star Formation Rate) [$M_\odot / yr$]",fontsize=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
plt.yscale('log')
plt.legend(loc='upper right')
plt.savefig("pls4.png")
plt.close()
'''
df1= pd.read_csv("tng33medfiltnorm.csv")
df2= pd.read_csv("tng67medfiltnorm.csv")
df3= pd.read_csv("tng99medfiltnorm.csv")
lower = 7
upper = 8
means = []
while upper <=11:
    df = pd.read_csv("tng99medfiltnorm.csv")
    df=df.dropna()
    df=df[df['mass']>lower]
    df=df[df['mass']<upper]
    slopes = list(df['slope2'])
    means.append(np.mean(slopes))
    print("for box {} to {}: mean slope is {}: N subhalos is {}".format(lower,upper, np.mean(slopes),len(slopes)))
    lower = upper
    upper= upper+1

print(means)

means33=    [-0.005246941276787649, -0.02190671465731498, -0.07324682074614279, -0.05564290989163309]
means67 =   [-0.005795415537161383, -0.019661913418380198, -0.06354791514996865, -0.05015920978972174]
means99 =   [-0.007283413953237877, -0.014826997854869945, -0.07107744409496307, -0.04343119472306738]

for i in range(4):
    diff1 = means99[i]-means67[i]
    diff2 = means67[i]-means33[i]
    totalchange = means99[i]-means33[i]
    print("diff1 = {}: diff2 = {}".format(diff1,diff2))
    print ("total change is {}".format(totalchange))
    print("pc total change = {}".format((totalchange/means33[i])*100))
    



plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

fig,axs = plt.subplots(nrows = 1, ncols = 3,figsize=(30,8),dpi=500)

divider0 = make_axes_locatable(axs[0])  ;cax0 = divider0.append_axes('right', size='5%', pad=0.05)
divider1 = make_axes_locatable(axs[1])  ;cax1 = divider1.append_axes('right', size='5%', pad=0.05)
divider2 = make_axes_locatable(axs[2])  ;cax2 = divider2.append_axes('right', size='5%', pad=0.05)

for i in range(3):
    axs[i].grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
    axs[i].tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
    axs[i].set_yscale('log')
    #axs[i].set_xscale('log')
    axs[i].set_xlabel("Total Mass [$M_\odot$]",fontsize=20)
    axs[i].set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
    

axs[0].set_title("Snapshot 33: z=2",fontsize=20)
axs[1].set_title("Snapshot 67: z=0.5",fontsize=20)
axs[2].set_title("Snapshot 99: z=0",fontsize=20)

im0=axs[0].scatter(df1['mass'],(df1['sfr']),c=(df1['slope2']),cmap='magma',vmin=-0.1,vmax=0.05)
im1=axs[1].scatter(df2['mass'],(df2['sfr']),c=(df2['slope2']),cmap='magma',vmin=-0.1,vmax=0.05)
im2=axs[2].scatter(df3['mass'],(df3['sfr']),c=(df3['slope2']),cmap='magma',vmin=-0.1,vmax=0.05)

fig.colorbar(im0, cax=cax0, orientation='vertical')
fig.colorbar(im1, cax=cax1, orientation='vertical')
fig.colorbar(im2, cax=cax2, orientation='vertical')

fig.tight_layout()
fig.subplots_adjust(top=0.89)
fig.suptitle("Redshift progression of galaxy Mass-SFR Relation", fontsize=20)
fig.savefig('inner_slopes_norm.png')

'''