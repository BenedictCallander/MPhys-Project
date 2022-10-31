import logging
from random import random # http logging for debugging purpouses
import time #runtime calculation import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import illustris_python as il


headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start = time.time()
#basePath='/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG100-1/output'
def get(path, params = None):
    #utility function for API reading 

    #Make API request - 
    # Path: url to api page 
    #Params - misc ; Headers = api key 
    r = requests.get(path, params=params, headers=headers)

    #HTTP code - raise error if code return is not 200 (success)
    r.raise_for_status()
    
    #detect content type (json or hdf5) - run appropriate download programme
    
    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

class visualisation:
    def __init__(self, df_g, df_s, subID, snapID, simID,crit_dist):
        self.df_g=pd.DataFrame({"x": 1, "y":2, "z":3,"rad":4,"m":10, "met":15})
        self.df_s=pd.DataFrame({"x": 1, "y":2, "z":3,"rad":4,"m":10, "met":15})


        self.df_g = df_g #gas dataframe
        self.df_s = df_s #stars dataframe
        self.subID = subID
        self.snapID = snapID
        self.simID = simID
        self.crit_dist = crit_dist
    def visual(self,type, quant, decp, annuli_pc):
        if (type=='gas'):
            df = self.df_g
            if(quant=='mass'):
                df_valid = df.round(decp)
                annul1= annuli_pc*self.crit_dist
                df_valid = self.df[self.df['rad']<annul1]
                df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
                plt.figure(figsize=(21,15))
                plt.style.use('dark_background')
                plt.scatter(-df_valid['x'],-df_valid['y'],c=(np.log10(df_valid['m'])),cmap='inferno', vmin=np.log10(min(df_valid['m'])), vmax = 0.95*np.log10(max(df_valid['m'])))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='log10(Gas Mass)')
                plt.title('Gas Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
                filename = 'temppng/Mgass_{}_sub_{}.png'.format(self.simID, self.subID)
                plt.savefig(filename)
                plt.close()
            elif(quant=='metallicity'):
                df_valid = df.round(decp)
                annul1= annuli_pc*self.crit_dist
                df_valid = self.df[self.df['rad']<annul1]
                df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
                plt.figure(figsize=(21,15))
                plt.style.use('dark_background')
                plt.scatter(-df_valid['x'],-df_valid['y'],c=(df_valid['met']),cmap='inferno', vmin=min(df_valid['m']), vmax = max(df_valid['m']))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='log10(Gas Metallicity)')
                plt.title('Metallicity Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
                filename = 'temppng/met_gas_{}_sub_{}.png'.format(self.simID, self.subID)
                plt.savefig(filename)
                plt.close()

        elif(type=='stars'):
            self.df = self.df_s
            if (quant=='mass'):
                df_valid = df.round(decp)
                annul1= annuli_pc*self.crit_dist
                df_valid = self.df[self.df['rad']<annul1]
                df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
                plt.figure(figsize=(21,15))
                plt.style.use('dark_background')
                plt.scatter(-df_valid['x'],-df_valid['y'],c=(np.log10(df_valid['m'])),cmap='inferno', vmin=np.log10(min(df_valid['m'])), vmax = 0.95*np.log10(max(df_valid['m'])))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='log10(Stellar Mass)')
                plt.title('Gas Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
                filename = 'temppng/Mstar_{}_sub_{}.png'.format(self.simID, self.subID)
                plt.savefig(filename)
                plt.close()
            
            elif(quant=='metallicity'):
                df_valid = df.round(decp)
                annul1= annuli_pc*self.crit_dist
                df_valid = self.df[self.df['rad']<annul1]
                df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
                plt.figure(figsize=(21,15))
                plt.style.use('dark_background')
                plt.scatter(-df_valid['x'],-df_valid['y'],c=(df_valid['met']),cmap='inferno', vmin=min(df_valid['m']), vmax = max(df_valid['m']))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='Metallicity)')
                plt.title('Metallicity Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
                filename = 'temppng/met_star_{}_sub_{}.png'.format(self.simID, self.subID)
                plt.savefig(filename)
                plt.close()

    def metgrad (self,type,decp, annuli_pc):
        if (type=='gas'):
            self.df = self.df_g
            df_valid = self.df.round(decp)
            annul1= annuli_pc*self.crit_dist
            df_valid = self.df[self.df['rad']<annul1]
            df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
            plt.figure(figsize=(21,15))
            plt.scatter(df_valid['rad'], (12+np.log10(df_valid['met'])), c=df_valid['m'], cmap = 'viridis')
            plt.xlabel('Radial Distance [kpc/h]')
            plt.ylabel('12+log10(O/H) [kpc/h]')
            plt.colorbar(label='Gass mass')
            plt.title('Gas Metallicity Gradient for SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
            filename = 'temppng/metgradGAS_{}_sub_{}.png'.format(self.simID, self.subID)
            plt.savefig(filename)
            plt.close()

        elif (type=='stars'):
            self.df = self.df_s
            df_valid = self.df.round(decp)
            annul1= annuli_pc*self.crit_dist
            df_valid = self.df[self.df['rad']<annul1]
            df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
            plt.figure(figsize=(21,15))
            plt.scatter(df_valid['rad'], (12+np.log10(df_valid['met'])), c=df_valid['m'], cmap = 'viridis')
            plt.xlabel('Radial Distance [kpc/h]')
            plt.ylabel('12+log10(O/H) [kpc/h]')
            plt.colorbar(label='Gass mass')
            plt.title('Stellar Metallicity gradient for SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
            filename = 'temppng/metgradSTAR_{}_sub_{}.png'.format(self.simID, self.subID)
            plt.savefig(filename)
            plt.close()


#"mygalaxy%d"%ID=Galaxy()