

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import module



final_tree_14 = np.loadtxt("tree_14_res.txt", dtype = np.array)
final_tree_13 = np.loadtxt("tree_13_res.txt")

print(final_tree_14)

fs = []
for j in range(4):
    for i in range(8):
        fs.append(np.array(module.pairing(final_tree_14[j][i], final_tree_13[j][i])[1::]))


print(type(fs))

a = np.array(fs)
print(type(a))
print(a)
rows_with_zero = np.any(a == 0, axis = 0)
print(rows_with_zero)
a = a[~rows_with_zero]
f = sum(np.array(fs).transpose())
f2 = sum(a.transpose())
k_i = np.array([0, 0, 1, 1, 1, 8, 0, 0, 0, 1, 0, 4, 0, 0, 0, 3, 0, 0, 0, 3])
mus = np.zeros
print(mus)
mus = sum(k_i)/f
mus2 = sum(k_i)/f2
print(mus)
#print(mus)
logL = []
#print('Valami',np.shape(a), a)
for i in range(15):
    logL.append(module.loglikelihood(mus2[i], np.array(a)[i], k_i )) 
print(logL)
L = new_list = [-15 if flag else logL.pop(0) for flag in rows_with_zero]
print(L)

# Define the values for x, y, and the 4th parameter
x_values = np.array([2, 4, 6, 8, 16, 32, 64])
y_values = np.array([1, 2, 3, 4])

# Generate the grid of x and y coordinates
X, Y = np.meshgrid(x_values, y_values)

# Define the z values (32 points for each x, y pair)
z_values = mus  # Example z values, replace with your actual data

# Define the 4th parameter values (32 points for each x, y pair)
parameter_values = L  # Example parameter values, replace with your actual data

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