import pandas as pd
import os

df = pd.read_csv("traceids.csv")
ids = list(df['id'])

for i in ids:
    dir_name = "files/binary/historycutouts/evdir_{}".format(i)
    os.makedirs(dir_name)
