#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 


df = pd.read_csv('bigdata1.csv')
print(df['sfr'])

plt.figure(figsize = (20,12))
plt.scatter(df['mass_stars'], 12+np.log(df['metallicity_stars']), marker = '+', c = df['id'], cmap = 'plasma')
plt.xscale('log')
plt.ylabel("metallicity_stars", fontsize = 25)
plt.xlabel("Mass)", fontsize = 25)
plt.tick_params(axis = 'both', which = 'both',direction = 'inout', length = 15, labelsize = 15)
plt.colorbar(label='arbitrary 3rd')
plt.savefig('sfr1.png')
plt.show()
plt.close()


plt.figure(figsize = (20,12))
plt.scatter(df['mass_gas'], 12+np.log(df['m_gas_sfr_weighted']), marker = '+',c = df['id'], cmap = 'hsv')
plt.xscale('log')
plt.ylabel("metallicity_stars", fontsize = 25)
plt.xlabel("Mass)", fontsize = 25)
plt.tick_params(axis = 'both', which = 'both',direction = 'inout', length = 15, labelsize = 15)
plt.colorbar(label='arbitrary 3rd')
plt.savefig('sfr2.png')
plt.show()
plt.close()



# %%
