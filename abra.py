import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

# Define the values for x, y, and the 4th parameter
x_values = np.array([2,6,10,14,18,22,26,30])
y_values = np.array([1, 2, 4, 8])

# Generate the grid of x and y coordinates
X, Y = np.meshgrid(x_values, y_values)

# Define the z values (32 points for each x, y pair)
z_values = np.array([0.31921068, 0.66747573, 0.71428571, 0.68365444, 0.79537238, 0.90163934,
 1.04562738, 1.17021277, 0.16267377, 0.4587156,  0.42868277, 0.46257359,
 0.52231719, 0.54563492, 0.59556037, 0.65243179, 0.07958327, 0.24509804,
 0.31187978, 0.30829596, 0.27281746, 0.31064671, 0.35132546, 0.36030134,
 0.04028714, 0.12066696, 0.20065669, 0.22231205, 0.21666338, 0.20109689,
 0.21938572, 0.24504344])  # Example z values, replace with your actual data

# Define the 4th parameter values (32 points for each x, y pair)
parameter_values = np.array([-15, -15, -4.322729916258087, -7.028497307769387, -8.001867074279442, -8.248189193142212, -11.29968347110716, -11.456386250067355, -15, -15, -15, -3.2715417418122676, -4.967931198544932, -5.872190240148164, -7.023820601921061, -11.523403312065424, -15, -15, -15, -15, -2.4750846230126697, -15, -4.494612122720774, -6.8345241636021115, -15, -15, -15, -15, -15, -15, -15, -5.802357721599233])  # Example 4th parameter values
  # Example parameter values, replace with your actual data

# Create a 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the scatter points
scatter = ax.scatter(X, Y, z_values, c=parameter_values, cmap='viridis')

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('3D Scatter Plot with Color')

# Create a ScalarMappable object for the color bar
sm = cm.ScalarMappable(cmap='viridis')
sm.set_array(parameter_values)
plt.colorbar(sm, label='4th Parameter')

# Show plot
plt.show()
