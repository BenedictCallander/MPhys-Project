import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def shapemeans(i):
    df = pd.read_csv("alldatadone.csv")
    df.dropna()
    df = df[df['snapshot'].isin([i])]
    slopes1 =list(df['slope1'])
    slopes2 =list(df['slope2'])
    slopes3 = list(df['slope3'])
    sfr = list(df['sfr'])
    print("mean slopes for snapshot {}".format(i))
    mean1 = np.mean(slopes1)
    mean2 = np.mean(slopes2)
    mean3 = np.mean(slopes3)
    mean4 = np.mean(sfr)
    print("slope1={} : slope2 ={} : slope3={}".format(mean1,mean2,mean3))
    print("mean SFR: {}".format(mean4))
    return(mean1,mean2,mean3,mean4)
snapshots=(21,33,50,67,78,91,99)

l1=[];l2=[];l3=[];l4=[]
for i in snapshots:
    mean1,mean2,mean3,mean4 = shapemeans(i)
    l1.append(mean1)
    l2.append(mean2)
    l3.append(mean3)
    l4.append(mean4)
    

plt.figure(figsize=(20,12))
plt.plot(snapshots,l1,'g-',label = 'Inner Gradient')
plt.plot(snapshots,l2,'r-',label = 'Centre Gradient')
plt.plot(snapshots,l3,'b-',label = 'Outer Gradient')
plt.plot(snapshots,l4,'k--',label='SFR')
plt.legend(loc='upper right')
plt.savefig("means.png")
plt.close()

    
    


