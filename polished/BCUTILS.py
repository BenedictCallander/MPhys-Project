#
# BCUTILS.py : utility containing function
#

#Plotting, numerical functions and dataset manipulation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pwlf

#Read data from web API and monitor HTTP traffic 
import requests  

#specialised functions to query Illustris TNG data 
import illustris_python as il

#Own module containing utility functions 
import BCUTILS

#runtime calculation 
import time

#Computational functions - simultaneous calculations to make use of multi-core CPU
from joblib import Parallel, delayed

#specific functions for fitting utilities
from scipy.optimize import curve_fit
from scipy.signal import medfilt, savgol_filter






class UTILITY:
    def line(m,x,b):
        y= m*x+b
        return y
    


class ANALYSIS:
    def get_signs(a, b, c):
        # Initialize an empty list
        signs = []
        
        # Append the signs of the input values to the list
        signs.append('negative' if a < 0 else 'positive' if a > 0 else 'zero')
        signs.append('negative' if b < 0 else 'positive' if b > 0 else 'zero')
        signs.append('negative' if c < 0 else 'positive' if c > 0 else 'zero')
        
        # Return the list of signs
        return signs
    
    def bootstrap_linear()