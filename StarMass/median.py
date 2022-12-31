import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
right_rad = np.linspace(0,10,101)

def getfname(i):
    fname = ("992slopes/slopedata_{}.csv".format(i))
    return fname
'''
dfin = pd.read_csv("tng99MS.csv")
dfin.mass = np.log10(dfin.mass)
ids = list(dfin['id'])


df_list = []
total_not_done = 0
for i in ids:
    try:
        fname = getfname(i)
        dftemp = pd.read_csv(fname)
        #dftemp = dftemp.insert(3,"radmed",right_rad,True)
        df_list.append(dftemp)
    except FileNotFoundError:
        print(i)
        total_not_done = total_not_done+1
print(total_not_done)
df = pd.concat(df_list)
median_plot = df.groupby('rad').mean()
median_plot.to_csv("med99.csv")
plt.plot(right_rad, median_plot['met'], 'r-')
plt.savefig("temp.png")
'''
df1 = pd.read_csv("med33.csv")
df2 = pd.read_csv("med67.csv")
df3 = pd.read_csv("med99.csv")
def line(m,x,b):
    y=(m*x)+b
    return y
right_rad = np.linspace(0,10,101)
popt,pcov = curve_fit(line, right_rad, (df1['met']))
popt1,pcov = curve_fit(line, right_rad, (df2['met']))
popt2,pcov = curve_fit(line, right_rad, (df3['met']))


df332 = pd.read_csv("high33.csv")
df33 = pd.read_csv("low33.csv")

df672 = pd.read_csv("high67.csv")
df67 = pd.read_csv("low67.csv")

df992 = pd.read_csv("high99.csv")
df99 = pd.read_csv("low99.csv")

fig,axs = plt.subplots(1,2,figsize=(20,8))
for i in range(2):
    axs[i].set_ylabel("12 + $log_{10}$(O/H)",fontsize=15)
    axs[i].set_xlabel("Radial Distance (Normalised Code Units)",fontsize=15)
    axs[i].grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
    axs[i].tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1,labelsize=15)
    
    
axs[0].plot(right_rad, (df33['met']),'r-',label = " $7<m<9.5  M_{\odot}$:z~2")
axs[0].plot(right_rad, (df332['met']),'y-',label = "$9.5<m<12  M_{\odot}$:z~2")

axs[0].plot(right_rad, (df67['met']),'b-',label = "$7<m<9.5  M_{\odot}$:z~0.5")
axs[0].plot(right_rad, (df672['met']),'g-',label = "$9.5<m<12  M_{\odot}$:z~0.5")


axs[0].plot(right_rad, (df99['met']),'k-',label = "$7<m<9.5  M_{\odot}$:z~0")
axs[0].plot(right_rad, (df992['met']),'c-',label = "$9.5<m<12  M_{\odot}$:z~0")

axs[0].set_title("Mass Binned Median Metallicity Profiles",fontsize=20)
axs[1].plot(right_rad,(df1['met']), 'r-',label="z~2: slope = {:.4f}".format(popt[0]))
axs[1].plot(right_rad,(df2['met']), 'b-',label="z~0.5 slope = {:.4f}".format(popt1[0]))
axs[1].plot(right_rad,(df3['met']), 'g-',label="z=0 slope = {:.4f}".format(popt2[0]))
axs[0].legend(loc="upper right")
axs[1].legend(loc="upper right")

axs[1].set_title("Median Metallicity Profiles for TNG50-1 Subhalo Population",fontsize=20)


fig.tight_layout()
fig.savefig("temp2.png")

