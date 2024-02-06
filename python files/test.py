import csv
import random
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
vec = [random.normalvariate(0,1),random.normalvariate(0,1),random.normalvariate(0,1)]
mag = np.dot(vec,vec)
S = [v/mag for v in vec]
def get_circle():
	L = np.hypot(S[0], np.hypot(S[1], S[2]))
	if S[0] > 0: sig = L
	else: sig = -L
	h = S[0]+sig
	beta = -1/(sig*h)
	f = beta*S[1]
	g = beta*S[2]
	return [f*h, 1+f*S[1], f*S[2]], [g*h, g*S[1], 1+g*S[2]]

fig = plt.figure()
ax = fig.add_subplot(projection ='3d')
ax.quiver(0,0,0, S[0], S[1], S[2])
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
U, V = get_circle()
theta_list = np.linspace(0,2*np.pi, 100)
for theta in theta_list: 
	pos = [(s+np.cos(theta)*u+np.sin(theta)*v) for s,u,v in zip(S,U,V)]
	ax.scatter(pos[0], pos[1], pos[2])

plt.show()