import csv
import numpy as np
import statsmodels
import matplotlib.pyplot as plt

raw = np.load('.\python files\tmp\Raw2023-07-05 13_16_32.454088.csv.npy')
I0 = raw[::50, 0]

plot_acf(I0, lags= 100)
plt.show()