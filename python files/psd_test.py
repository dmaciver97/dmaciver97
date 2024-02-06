import numpy as np
from scipy.signal import periodogram
import random
import matplotlib.pyplot as plt

dt = 1e-5
w_i = [0]
N = 100000
for i in range(N):
    w_i.append(random.normalvariate(0,1)+w_i[i])

step_sizes = np.logspace(0, 5, 10000)
msd_w = []
tau = [step*dt for step in step_sizes]
for step in step_sizes:
	n = int(N/step)
	di = int(step)
	tmp_w = []
	for i in range(0, n-1):
		tmp_w.append((w_i[i+di]-w_i[i])**2)
	msd_w.append(np.mean(tmp_w))

plt.plot(tau, msd_w)
plt.xscale('log')
plt.yscale('log')
plt.show()