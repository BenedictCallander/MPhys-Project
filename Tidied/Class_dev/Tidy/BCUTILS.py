import logging
from random import random # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il


plt.figure(figsize=(20,12), dpi=500)
plt.style.use('dark_background')
plt.scatter(-df_valid['x'],-df_valid['y'],c=(np.log10(df_valid['m'])),cmap='inferno', vmin=(min(np.log10(df_valid['m']))),vmax =(max(np.log10(df_valid['m']))))
plt.xlabel('$\Delta x$ [kpc/h]')
plt.ylabel('$\Delta y$ [kpc/h]')
plt.colorbar(label='log10(Gas Mass)')
.title('Gas Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
filename = 'temppng/Mgass_{}_sub_{}.png'.format(self.simID, self.subID)
plt.savefig(filename)
plt.close()