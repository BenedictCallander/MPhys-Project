import numpy as np 
import matplotlib.pyplot as plt 

x = [0,5,10]
y=[7,6,4]

fig = plt.figure(figsize=(8,5))
plt.plot(x,y,'r-')
plt.xlabel("X")
plt.ylabel("y")
plt.ylim(0,10)
plt.savefig("shape4.png")
plt.close()


print((7-6)/5)