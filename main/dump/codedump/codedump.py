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
