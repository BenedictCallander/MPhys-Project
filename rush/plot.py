import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt

'''
df = pd.read_csv("tng67.csv")
def line(m,x,b):
    y=(m*x) +b
    return y 
def line1(x):
    y=(1*x)-11
    return y 
df.mass = ((df.mass*10e10)/0.7)
def MSfilterup(dfin):
    df = dfin
    ids = list(df['id'])
    masses = list(np.log10(df['mass']))
    sfr = list(np.log10(df['sfr']))
    valids = []
    for i in range(len(ids)):
        value=line1((masses[i]))
        if value<((sfr[i])):
            valids.append(ids[i])
            print(i)
        else:
            continue
    return valids
valid1 = MSfilterup(df)
df2 = df[df['id'].isin(valid1)].copy()
print(df2)
df2.to_csv("tng67MS.csv")
xvals = np.linspace(6,13,100)
plt.figure(figsize=(20,15))
plt.plot(np.log10(df2['mass']),np.log10(df2['sfr']),'r+',zorder=1)
plt.plot(np.log10(df['mass']),np.log10(df['sfr']),'g+',zorder=0)
plt.plot(xvals,line1(xvals),'r--')
plt.xlim(9,13)
plt.savefig("MSFR67.png")
plt.close()
'''
df = pd.read_csv("tng67slopes.csv")

print(min(df['slope2']))
print(max(df['slope2']))
print(np.mean(df['slope2']))

plt.figure(figsize=(20,15))
plt.scatter(np.log10(df['mass']),np.log10(df['sfr']),c= df['slope2'], cmap = 'magma',vmin = -0.4, vmax=0.25)
plt.xlabel("Log(Total Mass) [$M_\odot$]",fontsize=20)
plt.ylabel("Log(Star Formation Rate) [$M_\odot / yr$]",fontsize=20)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
plt.xlim(8,12)
plt.colorbar().set_label(label="Metallicity Gradient (dex $Kpc^{-1}$)",size=15)
plt.tight_layout()
plt.savefig("slopes67.png")
plt.close()
