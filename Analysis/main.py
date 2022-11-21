import time  # runtime calculation import numpy as np #data handling
from random import random

import illustris_python as il
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests  # obtain data from API server
import seaborn as sns 
import BCUTILS
headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}
start =time.time()

dfin = pd.read_csv("csv/tng33subhalos.csv")
df2 = pd.read_csv("csv/tng33slopes.csv")

BCUTILS.MSfilter(dfin,df2,'csv/tng33MAIN.csv')
end = time.time()

print("runtime : {}".format(end-start))