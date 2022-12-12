import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

df = pd.read_csv("tng99MS2slopes.csv")
df=df.dropna()
print(min(df['slope2']))
print(max(df['slope2']))
print(np.mean(df['slope2']))

df.slope2 = (df.slope2-df.slope2.min())/(df.slope2.max()-df.slope2.min())

plt.figure(figsize=(20,12))
plt.scatter(df['mass'],df['sfr'],c=df['slope2'], cmap='magma',vmin=0.2,vmax=0.6)#,vmin=0.25,vmax=0.75)
plt.xlabel("Total Mass [$M_\odot$]",fontsize=20)
plt.yscale('log')
plt.ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)
plt.colorbar(label='Outer Slope')
plt.savefig("outerslopes.png")
plt.close()

print(df)