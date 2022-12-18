#%%
import logging
from random import random # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
from joblib import Parallel, delayed
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il

df = pd.read_csv("MZR99.csv")


plt.figure(figsize=(20,12))

plt.title("Metallicity Gradients of TNG50-1 Main Sequence Subhalos at z=0.5")
plt.xlabel("Total Mass [$M_\odot$]",fontsize=20)
plt.ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
plt.plot(df['mass'],12+np.log10(df['met']),'g+')
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.5,alpha =0.5)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig("slopechar67.png")
plt.show()
plt.close()
# %%
