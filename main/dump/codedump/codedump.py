#
# CODE DUMP -> LEAVE COMPLETED RUN STATEMENTS OR WORKING CODE TO BE SAVED 
#


'''
returns = Parallel(n_jobs= 20)(delayed(subhalo_slope_analysis)(i) for i in valid_id)
df3=pd.DataFrame(returns,columns=['linslope','brokenslope'])

linearslopes = list(df3['linslope'])
brokenslopes = list(df3['brokenslope'])

print("Linear min: {}   MAX: {}".format(min(linearslopes),max(linearslopes)))
print("Broken min: {}   MAX: {}".format(min(brokenslopes),max(brokenslopes)))


BCUTILS.MSfilter(dfin,df2,'csv/tng99MAIN.csv')
Key Numbers: 

'''

'''
returns = Parallel(n_jobs= 20)(delayed(subhalo_analysis)(i) for i in valid_id)
df2=pd.DataFrame(returns,columns=['slope','met','id','sfr','inside','outside'])
df2.insert(5,'mass', dfin['mass'],True)
df2.dropna()
df2.to_csv("csv/tng33MSslopes.csv")
'''
    
    
    
import os
import requests

# URL of the API that contains the file
api_url = "https://www.example.com/api/file"

# Name of the file
file_name = "my_file.txt"

# Path to the directory where the file should be saved
directory = "./data"

# Check if the file already exists in the directory
if not os.path.exists(os.path.join(directory, file_name)):
    # Request the file from the API
    response = requests.get(api_url)
    
    # Save the file to disk
    with open(os.path.join(directory, file_name), "wb") as f:
        f.write(response.content)
        
    print(f"Successfully downloaded {file_name} to {directory}.")
else:
    print(f"{file_name} already exists in {directory}.")

#
#
#
#

# List of downloaded filenames
filenames = ['file1.txt', 'file2.txt', 'file3.jpg']

# Create directories based on filenames
for filename in filenames:
    directory_name = filename.split('.')[0]  # Get the file name without the extension
    if not os.path.exists(directory_name):  # Check if the directory already exists
        os.makedirs(directory_name)  # Create the directory



class MyClass:
    def __init__(self, name):
        # Store the name of the object in an instance variable
        self.name = name

        # Use the __setattr__ method to dynamically create a new object
        # with the name stored in the `name` instance variable
        self.__setattr__(name, "hello world")


# Create an instance of MyClass with the name "my_object"
my_class = MyClass("my_object")

# Print the dynamically named object to verify that it has been created
print(my_class.my_object)  # Output: "hello world"



import matplotlib.pyplot as plt
import numpy as np

# Create some dummy data to plot
data1 = np.random.randn(100)
data2 = np.random.randn(100)
data3 = np.random.randn(100)

# Create a figure and axes for the plot
fig, ax = plt.subplots()

# Plot the first dataset
ax.plot(data1, label="Data 1")

# Set the title and labels for the plot
ax.set_title("Animated Background Data")
ax.set_xlabel("X axis")
ax.set_ylabel("Y axis")

# Add a legend to the plot
ax.legend()

# Function to update the background data and redraw the plot
def update(curr):
    # Check if the animation is at the last frame, and if so, stop it
    if curr == 100:
        a.event_source.stop()

    # Clear the current axes
    ax.clear()

    # Select the current dataset to plot based on the current frame
    if curr < 50:
        data = data1
    elif curr < 75:
        data = data2
    else:
        data = data3

    # Plot the selected dataset
    ax.plot(data, label="Data")

    # Set the title and labels for the plot
    ax.set_title("Animated Background Data")
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")

    # Add a legend to the plot
    ax.legend()

# Create an animation using the update function
a = animation.FuncAnimation(fig, update, interval=100)

# Show the plot
plt.show()

treepath = "files/binary/historycutouts/evdir_{}/treedata_{}.csv".format(self.primesub,self.primesub)
        df=pd.read_csv(treepath)




import numpy as np
from scipy.optimize import curve_fit

# Define the x-values for the piecewise function
x = np.array([0, 1, 2, 3, 4, 5])

# Define the y-values for the piecewise function
y = np.array([1, 3, 7, 9, 12, 15])

# Define the breakpoints for the piecewise function
# These are the initial guesses for the breakpoints
breakpoints = [1, 3]

# Define the piecewise function using the numpy.piecewise function
# This will define a piecewise linear function with two segments
# (one from x=0 to x=1, and another from x=3 to x=5)
pw_linear = np.piecewise(x, [x < breakpoints[0], x >= breakpoints[1]], [lambda x: x, lambda x: x - 2])

# Use the scipy.optimize.curve_fit function to find the optimal breakpoints
# that minimize the error between the piecewise linear function and the actual data
optimal_breakpoints, _ = curve_fit(pw_linear, x, y)

# Print the optimal breakpoints
print(optimal_breakpoints)



import numpy as np
import matplotlib.pyplot as plt

# Define the function to compute the piecewise linear fit with 0, 1, or 2 breakpoints
def piecewise_linear_fit(x, y, max_breakpoints=2):
    # Compute the initial least squares fit
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y)[0]

    # Compute the residuals from the initial fit
    y_fit = m * x + c
    residuals = y - y_fit

    # Check if we should add a breakpoint
    if max_breakpoints > 0:
        # Compute the maximum residual and the corresponding x-value
        max_residual = np.max(np.abs(residuals))
        x_break = x[np.argmax(np.abs(residuals))]

        # Compute the least squares fit with a breakpoint at the x_break value
        x1 = x[x < x_break]
        y1 = y[x < x_break]
        x2 = x[x >= x_break]
        y2 = y[x >= x_break]

        A1 = np.vstack([x1, np.ones(len(x1))]).T
        m1, c1 = np.linalg.lstsq(A1, y1)[0]
        A2 = np.vstack([x2, np.ones(len(x2))]).T
        m2, c2 = np.linalg.lstsq(A2, y2)[0]

        # Compute the residuals from the fit with a breakpoint
        y1_fit = m1 * x1 + c1
        y2_fit = m2 * x2 + c2
        y_break_fit = np.concatenate([y1_fit, y2_fit])
        residuals_break = y - y_break_fit

        # Check if the fit with a breakpoint is better
        if np.sum(residuals_break**2) < np.sum(residuals**2):
            # Use the fit with a breakpoint
            m = [m1, m2]
            c = [c1, c2]
            residuals = residuals_break

            # Check if we should add a second breakpoint
            if max_breakpoints > 1:
                # Compute the maximum residual and the corresponding x-value
                max_residual = np.max(np.abs(residuals))
                x_break = x[np.argmax(np.abs(residuals))]

                # Compute the least squares fit with a breakpoint at the x_break value
                x1 = x[x < x_break]
                y1 = y[x < x_break]
                x2 = x[x >= x_break]
                y2 = y[x >= x_break]

                A1 = np.vstack([x1, np.ones(len(x1))]).T
                
                



from numba import njit

@njit
def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

Yes, the code you provided could potentially be made more efficient and run faster. Here are a few suggestions to try:

    Use the @njit decorator from the numba library to compile the piecewise_linear() function to machine code, which can improve performance.

from numba import njit

@njit
def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

    Use the np.vectorize() function to vectorize the piecewise_linear() function, which can improve performance by applying the function to multiple input values simultaneously.

def piecewise_linear(x, x0, y0, k1, k2):
    return np.vectorize(lambda x: np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0]))(x)

    Use the np.isnan() function to check for NaN values in the x and y data before fitting the piecewise linear function, which can improve performance by avoiding unnecessary calculations.

def fit_piecewise_linear(x, y):
    # Check for NaN values in the data
    if np.isnan(x).any() or np.isnan(y).any():
        return

    # Define the initial values for the fitting parameters
    p0 = [0, 0, 0, 0]

    # Fit the piecewise linear function with 0 breakpoints
    res = optimize.curve_fit(piecewise_linear, x, y, p0)
    params = res[0]
    success = res[1]

    # If the fit with 0 breakpoints is not successful, try fitting with 1 breakpoint
    if success < 0.95:
        def piecewise_linear_1b(x, x1, y0, k1, k2):
            return np.piecewise(x, [x < x1], [lambda x:k1*x + y0-k1*x1, lambda x:k2*x + y0-k2*x1])
        p0 = [np.median(x), 0, 0, 0]
        res = optimize.curve_fit(piecewise_linear_1b, x, y, p0)
        params = res[0]
        success = res[1]

        # If the fit with 1 breakpoint is not successful, try fitting with 2 breakpoints
        if success < 0.95:
            def piecewise_linear_2b(x, x1, x2, y0, k1, k2, k3):
                return np.piecewise(x, [x < x1, (x >= x1) & (x < x2)], [lambda x:k1*x + y0-k1*x1, lambda x:k2*x + y0-k2*x2, lambda x:k3*x +


import numpy as np
from scipy import spatial

# Define a function to align the subhalo to be perpendicular with the z direction
def align_subhalo(pos, vel):
    # Calculate the center of mass of the subhalo
    com = np.mean(pos, axis=0)

    # Translate the positions and velocities to the center of mass frame
    pos = pos - com
    vel = vel - np.mean(vel, axis=0)

    # Calculate the principal axes of the subhalo using the inertia tensor
    I = np.sum(pos**2, axis=1)
    inertia_tensor = np.array([[np.sum(pos[:,0]**2), np.sum(pos[:,0]*pos[:,1])],
                               [np.sum(pos[:,0]*pos[:,1]), np.sum(pos[:,1]**2)]])
    eig_vals, eig_vecs = np.linalg.eig(inertia_tensor)

    # Align the subhalo to be perpendicular with the z direction
    pos = np.dot(pos, eig_vecs)
    vel = np.dot(vel, eig_vecs)

    # Return the aligned positions and velocities
    return pos, vel

# Define a function to compute the metallicity gradient
def metallicity_gradient(pos, met):
    # Compute the distance of each particle from the center of the subhalo
    dist = np.sqrt(pos[:,0]**2 + pos[:,1]**2)

    # Compute the median metallicity of the subhalo as a function of distance
    median_met = []
    for d in np.unique(dist):
        median_met.append(np.median(met[dist == d]))

    # Fit a linear function to the median metallicity as a function of distance
    gradient, intercept = np.polyfit(dist, median_met, deg=1)

    # Return the metallicity gradient
    return gradient

# Load the particle data for the subhalo
pos = np.load('positions.npy')
vel = np.load('velocities.npy')
met = np.load('metallicities.npy')

# Align the subhalo to be perpendicular with the z direction
pos, vel = align_subhalo(pos, vel)

# Compute the metallicity gradient
gradient = metallicity_gradient(pos, met)
print(gradient)


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


import pandas as pd
import numpy as np

# Define a function to compute the square root of the sum of the squares of two columns
def compute_sqrt_sum_squares(row, col1, col2):
    return np.sqrt(row[col1]**2 + row[col2]**2)

# Load the data into a DataFrame
df = pd.read_csv('data.csv')

# Compute the square root of the sum of the squares of two columns and add the result as a new column
df['sqrt_sum_squares'] = df.apply(compute_sqrt_sum_squares, args=('col1', 'col2'), axis=1)


import pandas as pd

# Load the DataFrame
df = pd.read_csv('data.csv')

# Define the limit
limit = 10

# Use the apply() method to filter the rows
filtered_df = df.apply(lambda x: x['z'] < limit, axis=1)

# View the filtered DataFrame
print(filtered_df)



import os

# Define the base directory
base_dir = '/path/to/base/directory'

# Iterate over the range of numbered directories
for i in range(1, 101):
    # Construct the path to the directory
    dir_path = os.path.join(base_dir, f'dir_{i:03d}')

    # Check if the directory exists
    if os.path.exists(dir_path):
        # Iterate over the 6 numbered files in the directory
        for j in range(1, 7):
            # Construct the path to the file
            file_path = os.path.join(dir_path, f'file_{j:03d}.txt')

            # Check if the file exists
            if os.path.exists(file_path):
                # Process the file
                with open(file_path, 'r') as f:
                    # Read the contents of the file
                    contents = f.read()

                    # Do something with the contents of the file
                    ...



import numpy as np
from scipy import optimize

def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def fit_piecewise_linear(x, y):
    # Define the initial values for the fitting parameters
    p0 = [0, 0, 0, 0]

    # Fit the piecewise linear function with 0 breakpoints
    res = optimize.curve_fit(piecewise_linear, x, y, p0)
    params = res[0]
    success = res[1]

    # If the fit with 0 breakpoints is not successful, try fitting with 1 breakpoint
    if success < 0.95:
        def piecewise_linear_1b(x, x1, y0, k1, k2):
            return np.piecewise(x, [x < x1], [lambda x:k1*x + y0-k1*x1, lambda x:k2*x + y0-k2*x1])
        p0 = [np.median(x), 0, 0, 0]
        res = optimize.curve_fit(piecewise_linear_1b, x, y, p0)
        params = res[0]
        success = res[1]

        # If the fit with 1 breakpoint is not successful, try fitting with 2 breakpoints
        if success < 0.95:
            def piecewise_linear_2b(x, x1, x2, y0, k1, k2, k3):
                return np.piecewise(x, [x < x1, (x >= x1) & (x < x2)], [lambda x:k1*x + y0-k1*x1, lambda x:k2*x + y0-k2*x2, lambda x:k3*x + y0-k3*x2])
            p0 = [np.min(x), np.max(x), 0, 0, 0, 0]
            res = optimize.curve_fit(piecewise_linear_2b, x, y, p0)
            params = res[0]

    # Return the fitting parameters
    return params

# Define the x and y data
x = np.array([1, 2, 3, 4, 5, 6])
y = np.array([1, 2, 1.3, 3, 2, 5])

# Fit the piecewise linear function to the data
params = fit_piecewise_linear(x, y)
print(params)



import numpy as np
import matplotlib.pyplot as plt

# Define the function to compute the broken linear fit
def broken_linear_fit(x, y, breakpoint):
    # Compute the least squares fit of the data to a linear function before the breakpoint
    x1 = x[x < breakpoint]
    y1 = y[x < breakpoint]
    A1 = np.vstack([x1, np.ones(len(x1))]).T
    m1, c1 = np.linalg.lstsq(A1, y1)[0]

    # Compute the least squares fit of the data to a linear function after the breakpoint
    x2 = x[x >= breakpoint]
    y2 = y[x >= breakpoint]
    A2 = np.vstack([x2, np.ones(len(x2))]).T
    m2, c2 = np.linalg.lstsq(A2, y2)[0]

    # Compute the fitted values at the breakpoint
    y1_fit = m1 * breakpoint + c1
    y2_fit = m2 * breakpoint + c2

    # Compute the residuals at the breakpoint
    residual = y1_fit - y2_fit

    return m1, c1, m2, c2, residual

# Load the data into a NumPy array
data = np.loadtxt('data.txt')

# Extract the x and y values from the array
x = data[:, 0]
y = data[:, 1]

# Set the range of breakpoints to consider
breakpoints = np.arange(x[0], x[-1], 0.01)

# Compute the broken linear fits for each breakpoint
fits = [broken_linear_fit(x, y, bp) for bp in breakpoints]

# Find the breakpoint with the

