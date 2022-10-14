import logging # http logging for debugging purpouses
import time #runtime calculation 
import numpy as np #data handling 
import requests #obtain data from API server
import h5py #binary file manipulation
import pandas as pd 
import matplotlib.pyplot as plt 
import illustris_python as il
basepath = '/x/Physics/AstroPhysics/Shared-New/DATA/IllustrisTNG/TNG50-1/output'
fields = il.groupcat.loadsubhalos