import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#^ define particle path and parameters
radius_p = 1.57e-6#float(input("Particle radius [microns]: "))*1e-6
path_radius = 4.71e-6#float(input("Orbit radius [microns]: "))*1e-6
omega_0 = 1#float(input("Orbital frequency [Hz]: "))
period = 1/omega_0

#^ compute velocity 'U' along path
fig, ax  = plt.subplots()
scat, = ax.plot([],[], 'ro')
vec = ax.arrow([],[], [], [])

def init():
    ax.set_xlim(-path_radius-radius_p, path_radius+radius_p)
    ax.set_ylim(-path_radius-radius_p, path_radius+radius_p)
    return scat,

def update(frame):
    x_data=path_radius*np.cos(2*np.pi*frame)
    y_data=path_radius*np.sin(2*np.pi*frame)
    vx_data =-path_radius*2*np.pi/period*np.sin(2*np.pi*frame)
    vy_data =path_radius*2*np.pi/period*np.cos(2*np.pi*frame) 
    scat.set_data(x_data, y_data)
    vec.set_data(x = x_data, y = y_data, dx = vx_data*1e-2, dy= vy_data*1e-2,
                 width = 1e-8, head_width=0.1e-6, head_length=0.1e-6)

    return scat, vec

ani = FuncAnimation(fig, update, frames=np.linspace(0, 1, 100),
                    init_func=init, blit=True, interval = 10)

plt.show()