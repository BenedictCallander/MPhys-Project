#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 


df = pd.read_csv('bigdata2.csv')
print(df['sfr'])

plt.figure(figsize = (20,12))
plt.scatter(df['mass_stars'], 12+np.log(df['metallicity_stars']), marker = '+', c = np.log(df['sfr']), cmap = 'plasma')
plt.xscale('log')
plt.yscale('log')
plt.ylabel("Metallicity(Stars)", fontsize = 25)
plt.xlabel("Mass (Stars)", fontsize = 25)
plt.tick_params(axis = 'both', which = 'both',direction = 'inout', length = 15, labelsize = 15)
plt.colorbar(label='SFR')
plt.savefig('sfr1.png')
plt.show()
plt.close()


plt.figure(figsize = (20,12))
plt.scatter(df['mass_gas'], 12+np.log(df['metallicity_gas']), marker = '+',c = np.log(df['sfr']), cmap = 'hsv')
plt.xscale('log')
plt.ylabel("Metallicity(Gas)", fontsize = 25)
plt.xlabel("Mass(Gas))", fontsize = 25)
plt.tick_params(axis = 'both', which = 'both',direction = 'inout', length = 15, labelsize = 15)
plt.colorbar(label='arbitrary 3rd')
plt.savefig('sfr2.png')
plt.show()
plt.close()



# %%
