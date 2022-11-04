
# Simple Pendulum Simulation
# Benedict Callander : sn 199194504
#

#Import Libraries

import numpy as np 
from scipy import integrate 
import matplotlib.pyplot as plt

# define constants and values 

m = 1 #mass = 1kg
L = 1 #Length
b = float(input(" enter damping constant : ")) #damping constant value
g = 9.81 # value for g 
d_t = 0.02 # time step size for simulation
t_max = 20 # simulation duration

# conditions

th1_0 = np.pi/2 # initial angle
th2_0 = 0 # initial angular velocity
th_init=(th1_0, th2_0)

#create time steps 

t= np.arange(0,t_max, d_t)

def pendulum_sim(th_init, t, L=1, m=1, g=9.81):
    theta_dot_1 = th_init[1]
    theta_dot_2 = -b/m*th_init[1] - g/L*np.sin(th_init[0])
    return theta_dot_1, theta_dot_2

theta_vals_int = integrate.odeint(pendulum_sim, th_init, t)

tdotvalsp = (theta_vals_int[:,1]**2)/2
tvalsp = -(g/L) * np.cos(theta_vals_int[:,0]) 

tdotvals = theta_vals_int[:,1]
tvals = theta_vals_int[:,0]
plt.figure(figsize =(10,10))
plt.plot(tdotvals,tvals, 'r-')
plt.tick_params(direction='in',      
                length=7,            
                bottom='on',         
                left='on',
                top='on',
                right='on',
                
               )
plt.rcParams.update({'font.size':20})
plt.xlabel("x")
plt.ylabel("y")
plt.show()
plt.close()

#print(theta_vals_int)
plt.figure(figsize=(15,10))
plt.figsize=(50,50)
plt.plot(t, theta_vals_int[:,1],  "r-", label='Angular Velocity ('r'$\theta^.)$')
plt.plot(t, theta_vals_int[:,0],  "b-", label='Magnitude ('r'$\theta)$')
plt.legend()
plt.xlabel("Time(s)")
plt.ylabel('Magnitude ('r'$\theta)$')
plt.savefig("pendulum.png")