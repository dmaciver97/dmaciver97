import csv
import random
from tkinter import X
import numpy as np
import matplotlib.pyplot as plt

x_list = np.linspace(-1,1,10)
y_list = np.linspace(-1,1,10)

xx, yy = np.meshgrid(x_list, y_list)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)

for x, y in zip(xx, yy):
    phi = np.arctan2(x,y)
    r = (x**2+y**2)**0.5
    theta = np.arcsin(r)
    s = [r*np.sin(theta)*np.cos(phi),
         r*np.sin(theta)*np.sin(phi),
         r*np.cos(theta)]
    ax.scatter(s[0], s[1], s[2])

plt.show()