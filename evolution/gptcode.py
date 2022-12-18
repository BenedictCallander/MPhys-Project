import os
import h5py

# Define the lists
list1 = [1, 2, 3, 4]
list2 = ['a', 'b', 'c', 'd']

# Define the directory to save the file to
directory = '/path/to/directory'

# Create the directory if it does not already exist
if not os.path.exists(directory):
    os.makedirs(directory)

# Open the HDF5 file in write mode
with h5py.File(os.path.join(directory, 'lists.hdf5'), 'w') as f:
    # Write the lists to the HDF5 file
    f['list1'] = list1
    f['list2'] = list2

    # Set the keywords for the lists
    f['list1'].attrs['keyword'] = 'Numbers'
    f['list2'].attrs['keyword'] = 'Letters'
'''
To apply a 2-breakpoint piecewise linear fit to a dataset in Python, you can use the scipy.optimize.curve_fit
function from the SciPy library. This function allows you to fit a piecewise linear function to a set of data points.

Here is an example of how you can use the scipy.optimize.curve_fit function to apply a 2-breakpoint piecewise linear fit to a dataset:

'''

import numpy as np
from scipy.optimize import curve_fit

# Define the x and y data for the dataset
xdata = [0, 1, 2, 3, 4, 5]
ydata = [1, 3, 2, 5, 4, 6]

# Define a function that returns a piecewise linear function with 2 breakpoints
def piecewise_linear(x, x1, y1, x2, y2):
    return np.piecewise(x, [x < x1, (x >= x1) & (x < x2), x >= x2], [lambda x: x*y1, lambda x: y1 + (x-x1)*(y2-y1)/(x2-x1), lambda x: y2 + (x-x2)*(y2-y1)/(x2-x1)])

# Use the curve_fit function to fit a piecewise linear function to the data
params, covariance = curve_fit(piecewise_linear, xdata, ydata)

# Print the fitted parameters
print(params)

'''
This code will fit a 2-breakpoint piecewise linear function to the xdata and ydata provided, and it will print
the fitted parameters (i.e. the breakpoint coordinates) to the console. You can then use these parameters to plot the piecewise
linear function and visualize the fit.

Note that this code uses the numpy.piecewise function to define the piecewise linear function. This function allows you to specify different linear
functions for different ranges of the input x. In this case, the piecewise linear function has 3 pieces, with the first and third pieces being defined
by the slope and y-intercept of the overall line, and the second piece being defined by the slope and y-intercept between the two breakpoints.
You can adjust this function to define a piecewise linear function with a different number of breakpoints, or with different behavior for each piece.
'''


'''
Yes, the code to apply a 2-breakpoint piecewise linear fit to a dataset can be made faster for large datasets.
One way to do this is to use the scipy.optimize.curve_fit function with the method parameter set to lm. 
This will use the Levenberg-Marquardt optimization algorithm to fit the piecewise linear function to the data,
which is generally faster than the default optimization algorithm used by curve_fit.

Here is an example of how you can use the scipy.optimize.curve_fit function with the method parameter set to lm to apply
a 2-breakpoint piecewise linear fit to a large dataset:
'''


import numpy as np
from scipy.optimize import curve_fit

# Define the x and y data for the dataset
xdata = [0, 1, 2, 3, 4, 5]
ydata = [1, 3, 2, 5, 4, 6]

# Define a function that returns a piecewise linear function with 2 breakpoints
def piecewise_linear(x, x1, y1, x2, y2):
    return np.piecewise(x, [x < x1, (x >= x1) & (x < x2), x >= x2], [lambda x: x*y1, lambda x: y1 + (x-x1)*(y2-y1)/(x2-x1), lambda x: y2 + (x-x2)*(y2-y1)/(x2-x1)])

# Use the curve_fit function with the lm method to fit a piecewise linear function to the data
params, covariance = curve_fit(piecewise_linear, xdata, ydata, method="lm")

# Print the fitted parameters
print(params)

'''
This code will use the Levenberg-Marquardt optimization algorithm to fit the 2-breakpoint piecewise linear function to the xdata and ydata provided,
and it will print the fitted parameters to the console. You can then use these parameters to plot the piecewise linear function and visualize the fit.

Note that the Levenberg-Marquardt optimization algorithm may not always provide the best fit to the data,
so it is recommended to try using the default optimization algorithm as well to compare the results and see which one is more accurate
for your particular dataset.
'''


import numpy as np
from scipy.optimize import curve_fit

# Define the x and y data for the dataset
xdata = [0, 1, 2, 3, 4, 5]
ydata = [1, 3, 2, 5, 4, 6]

# Define a function that returns a piecewise linear function with 2 breakpoints, and also returns the slopes of each piece
def piecewise_linear(x, x1, y1, x2, y2):
    # Calculate the slopes of each piece
    m1 = y1 / x1
    m2 = (y2 - y1) / (x2 - x1)
    m3 = (y2 - y1) / (x2 - x1)
    
    # Return the piecewise linear function and the slopes
    return np.piecewise(x, [x < x1, (x >= x1) & (x < x2), x >= x2], [lambda x: x*m1, lambda x: y1 + (x-x1)*m2, lambda x: y2 + (x-x2)*m3]), m1, m2, m3

# Use the curve_fit function to fit a piecewise linear function to the data
params, covariance = curve_fit(piecewise_linear, xdata, ydata)

# Print the fitted parameters and the slopes
print(params)
print(piecewise_linear(xdata, *params)[1:])
'''
This code will fit a 2-breakpoint piecewise linear function to the xdata and ydata provided,
and it will print the fitted parameters and the slopes of each piece to the console.
You can then use these parameters and slopes to analyze the behavior of the piecewise linear function and visualize the fit.

Note that in this modified piecewise_linear function, the slopes of each piece are calculated using the fitted parameters
x1, y1, x2, and y2. These slopes are then returned as additional output parameters, along with the piecewise linear function itself.
You can adjust this code to return the slopes in a different format, or to calculate the slopes in a different way, depending on your specific needs.
'''


import numpy as np

# Sample scattered data
x = np.array([0, 1, 2, 3, 4, 5])
y = np.array([0, 1, 4, 9, 16, 25])

# Weights for each datapoint
weights = np.array([1, 2, 3, 2, 1, 2])

# Fit linear regression
A = np.vstack([x, np.ones(len(x))]).T
m, c = np.linalg.lstsq(A * weights[:, np.newaxis], y * weights, rcond=None)[0]

# Print result
print(f"y = {m} * x + {c}")



import matplotlib.pyplot as plt
import numpy as np

# Set the radius of the annulus
radius = 5

# Generate data for the x and y axes of the plot
theta = np.linspace(0, 2*np.pi, 100)
x = radius*np.cos(theta)
y = radius*np.sin(theta)

# Create a figure and axes object for the plot
fig, ax = plt.subplots()

# Plot the annulus using the `fill_between` function
ax.fill_between(x, y, radius*2, alpha=0.5)

# Show the plot
plt.show()









df_valid = df.round(decp)
                annul1= annuli_pc*self.crit_dist
                df_valid = df[df['rad']<annul1]
                df_valid = df_valid.groupby(['x','y'])['m'].sum().reset_index()
                plt.figure(figsize=(20,12), dpi=500)
                plt.style.use('dark_background')
                plt.scatter(-df_valid['x'],-df_valid['y'],c=(np.log10(df_valid['m'])),cmap='inferno', vmin=(min(np.log10(df_valid['m']))),vmax =(max(np.log10(df_valid['m']))))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='log10(Gas Mass)')
                plt.title('Gas Density of SubID {}: {} snapshot {}'.format(self.subID, self.simID, self.snapID))
                filename = 'suppng/Mgass_{}_sub_{}.png'.format(self.simID, self.subID)
                plt.savefig(filename)
                plt.close()
                
                
        
import numba

@numba.jit
def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

# use the piecewise_linear function as before
