import numpy as np
import matplotlib.pyplot as plt
import quaternion
from quaternion import as_spherical_coords

with open('C:\\Users\\xbb16146\\dmaciver97\\python_files\\sphere_in_LCP_Quadrants.txt', 'r') as in_file:
    lines = in_file.readlines()
    x_dat, y_dat, z_dat, = [], [], []
    q1, q2, q3, q4 = [], [], [], []
    index = []
    for line in lines[1:-1]:
        items = line.split(',')
        x_dat.append(float(items[0])+ float(items[1])- float(items[2]) - float(items[3])) 
        y_dat.append(float(items[0])+ float(items[2])- float(items[1]) - float(items[3]))
        z_dat.append(float(items[0])+ float(items[1])+ float(items[2]) + float(items[3]))
        q1.append(float(items[0]))
        q2.append(float(items[1]))
        q3.append(float(items[2]))
        q4.append(float(items[3]))
        index.append(int(float(items[4])/1e-5))

traj = np.load('Sphere_in_LCP_Trajectory.npy')
x_disp, y_disp, z_disp = [], [], []
for i in index:
    x_disp.append(traj[i,0])
    y_disp.append(traj[i,1])
    z_disp.append(traj[i,2])
    q = np.quaternion(traj[i,3],traj[i,4],traj[i,5],traj[i,6])
    
plt.scatter(x_disp, x_dat)
plt.show()



