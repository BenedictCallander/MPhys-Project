import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

df = pd.read_csv("csv/tng99MSslopes.csv")
df=df.dropna()
print(min(df['slope3']))
print(max(df['slope3']))
print(np.mean(df['slope3']))
df=df[df['slope3']>-3]
plt.figure(figsize=(20,12))
plt.scatter(df['mass'],df['sfr'],c=df['slope3'], cmap='magma',vmin=-20,vmax=1)
plt.xlabel("Total Mass [$M_\odot$]",fontsize=20)
plt.yscale('log')
plt.ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
plt.colorbar(label='Outer Slope')
plt.savefig("outerslopes.png")
plt.close()

print(df)