import logging # http logging for debugging purpouses
import time #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 
 
 
galaxy_df = pd.read_csv("test1.csv")

valid_galaxies = galaxy_df[galaxy_df['mass']<9.5]
print(len(valid_galaxies['id']))
