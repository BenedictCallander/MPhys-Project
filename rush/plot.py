import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt

df = pd.read_csv("csv/broken99.csv")
'''
plt.figure(figsize=(20,15))
plt.plot(np.log10(df.mass), df.slope, 'r+')
plt.plot(np.log10(df.mass), (df.slope+df.change),'g+')
plt.ylim(-0.03,0.03)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
#plt.xlabel("Log(Total Mass) [$M_\odot$]",fontsize=20)
#plt.ylabel("Log(Star Formation Rate) [$M_\odot / yr$]",fontsize=20)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1,labelsize=15)
plt.savefig("temp.png")
plt.close()
'''

#print(np.mean(df['change']))
print(min(df['slope2']))
print(max(df['slope2']))
print(np.mean(df['slope2']))

plt.figure(figsize=(20,15))
plt.scatter(np.log10(df['mass']),np.log10(df['sfr']),c= df['slope2'], cmap = 'magma',vmin = -0.08, vmax=0.05)
#plt.scatter(np.log10(df['mass']),np.log10(df['sfr']),c= df['slope2'], cmap = 'magma',vmin = -0.1, vmax=0.08)
plt.xlabel("Log(Total Mass) [$M_\odot$]",fontsize=20)
plt.ylabel("Log(Star Formation Rate) [$M_\odot / yr$]",fontsize=20)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
#plt.xlim(8.5,11.5)
#plt.ylim(-2,1)
plt.tight_layout()
plt.colorbar().set_label(label="Metallicity Gradient (dex $Kpc^{-1}$)",size=15)
plt.savefig("csv/ffs_99.png")
plt.close()

'''
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
df2.to_csv("tng33MS.csv")
xvals = np.linspace(6,13,100)
plt.figure(figsize=(20,15))
plt.plot(np.log10(df2['mass']),np.log10(df2['sfr']),'r+',zorder=1)
plt.plot(np.log10(df['mass']),np.log10(df['sfr']),'g+',zorder=0)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
plt.xlabel("Log(Total Mass) [$M_\odot$]",fontsize=20)
plt.ylabel("Log(Star Formation Rate) [$M_\odot / yr$]",fontsize=20)
plt.plot(xvals,line1(xvals),'r--')
plt.xlim(9,13)
plt.savefig("MSFR33.png")
plt.close()
'''

'''
#df = pd.read_csv("csv/tng99broken.csv")
df = pd.read_csv("csv/tng33bkpc.csv")
#df.slope2 =-1*abs(df.slope2)

#print(min(df['slope2']))
#print(max(df['slope2']))
#print(np.mean(df['slope2']))

print(min(df['slope2']))
print(max(df['slope2']))
print(np.mean(df['slope2']))

#print(df2)
#print(df)
plt.figure(figsize=(20,15))
plt.scatter(np.log10(df['mass']),np.log10(df['sfr']),c= df['slope2'], cmap = 'magma',vmin = -0.1, vmax=0.05)
#plt.scatter(np.log10(df['mass']),np.log10(df['sfr']),c= df['slope2'], cmap = 'magma',vmin = -0.1, vmax=0.08)
plt.xlabel("Log(Total Mass) [$M_\odot$]",fontsize=20)
plt.ylabel("Log(Star Formation Rate) [$M_\odot / yr$]",fontsize=20)
plt.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(visible=True,which='both',axis='both',color='grey',linestyle='-',linewidth=0.7,alpha =1)
plt.xlim(8.5,11.5)
plt.ylim(-2,1)
plt.tight_layout()
plt.colorbar().set_label(label="Metallicity Gradient (dex $Kpc^{-1}$)",size=15)
plt.savefig("csv/bkpc33.png")
plt.close()


df1 = pd.read_csv("csv/tng33bkpc.csv")
df2 = pd.read_csv("csv/tng67bkpc.csv")
df3 = pd.read_csv("csv/tng99bkpc.csv")


print(np.mean(df1['slope2']))
print(np.mean(df2['slope2']))
print(np.mean(df3['slope2']))
print(min(df1['slope2']))
print(min(df2['slope2']))
print(min(df3['slope2']))

print(max(df1['slope2']))
print(max(df2['slope2']))
print(max(df3['slope2']))
'''