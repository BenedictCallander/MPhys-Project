#%%
import logging # http logging for debugging purpouses
import time
from turtle import color #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 


#x = np.arange(-213,214,1)
#y = np.arange(-213,214,1)
#data=pd.read_csv('filet.csv')
data=pd.read_csv('/home/AstroPhysics-Shared/PROJECTS/MPhys_Schady/Projects/TNGmetgrads/BC763/filet.csv')
data.sort_values(['x','y','z'], inplace=True)
data2 = data.groupby(['x','y'])['m'].sum().reset_index()
#data2.to_csv('run.csv')
#print(data2)
critical_radius = np.arange(0,213,10)
zeros = np.zeros(len(critical_radius))

plt.figure(figsize=(21,15))
plt.style.use('dark_background')
#plt.hist2d(data2['x'],data2['y'], weights=data2['m'],bins=[500,500],cmap = 'inferno',vmin=(min(data2['m'])),vmax = max(data2['m']))
plt.scatter(data2['x'],data2['y'],c=np.log10(data2['m']),cmap='inferno',vmin=1.01*(min(np.log10(data2['m']))), vmax=0.99*(max(np.log10(data2['m']))))
#plt.plot(critical_radius,zeros,'g+', markersize=10,label='critical radius')
plt.xlabel('$\Delta x$ [kpc/h]')
plt.ylabel('$\Delta y$ [kpc/h]')
#plt.xlim(-100,100)
#plt.ylim(-100,100)
plt.colorbar(label='log10(Gas Density)')
plt.title('Gas Density of SubID 15129: TNG100-1 snapshot 70')
#plt.legend(loc='upper right')
plt.savefig('group15129.png')

plt.close()


# %%
