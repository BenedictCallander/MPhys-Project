import pandas as pd 
import os 


idslist = pd.read_csv("remove.csv")
ids=list(idslist['ids'])
filenames = []
for i in ids:
    filename = ("radpng/lin_fit_{}.png".format(i))
    filenames.append(filename)
    
for j in filenames:
    try:
        os.remove(j)
        print('removed {}'.format(j))
    except OSError:
        print('Delete failed')

