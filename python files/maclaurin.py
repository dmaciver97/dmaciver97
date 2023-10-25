import numpy as np
import matplotlib.pyplot as plt

x_range = np.linspace(-5, 5, 50)


for x in x_range:
    macsum = 1
    kterm = 1
    k = 1
    while kterm > 10e-11:
        kterm *= (x/k)
        macsum += kterm
        k += 1
    plt.scatter(x, macsum, c='black', s =3)


plt.plot(x_range, [np.exp(x) for x in x_range])
plt.vlines(0, -10, np.exp(5), color='black', linestyle ='--')
plt.show()
plt.close()
