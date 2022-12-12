import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

df = pd.read_csv("try672.csv")

plt.figure(figsize=(20,12))
plt.scatter(df['mass'],df['sfr'],c=df['Values'], cmap='tab20')#,vmin=0.2,vmax=0.6)#,vmin=0.25,vmax=0.75)
plt.xlabel("Total Mass [$M_\odot$]",fontsize=20)
plt.yscale('log')
plt.ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
plt.colorbar(label='Outer Slope')
plt.savefig("try672.png")
plt.close()

print(df)


plt.figure(figsize=(20,12))
plt.hist2d(df['Values'],df['mass'],bins=[100,100], weights=-df['sfr'],cmap ='tab20')
plt.xlabel("Shape Character")
plt.yscale('log')
plt.ylabel("Total Mass [$M_\odot$]",fontsize=20)
plt.savefig("see2.png")
plt.close()
