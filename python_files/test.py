import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import periodogram

with open('C:\\Users\\xbb16146\\dmaciver97\\python_files\\dimer_files\\data_07-May-2024 1to2QPD.txt', 'r') as infile:
    lines = infile.readlines()
    data = {'Time': [], 'X-dat': [], 'Y-dat': [], 'Z-dat': []}
    for line in lines[1:-1]:
        row = line.split(',')
        data['Time'].append(float(row[3][:-1]))
        data['X-dat'].append(float(row[0]))
        data['Y-dat'].append(float(row[1]))
        data['Z-dat'].append(float(row[2]))
    df = pd.DataFrame(data)
    pd.to_pickle(df, 'test.pkl')

plt.plot(df['Time'], df['X-dat'], label = r'$QPD_x$')
plt.plot(df['Time'], df['Y-dat'], label = r'$QPD_y$')
plt.xlabel('Time [s]')
plt.legend()
plt.ylabel('QPD Signal [V]')
plt.show()

exit(0)

kT = 1.38e-23*298
dt = df['Time'][1]-df['Time'][0]

def potential_well(data):    
    x_eq = np.mean(data)
    def p_x(x, kappa):
        return np.sqrt(kappa/(2*np.pi*kT))*np.exp(-(kappa/(2*kT))*(x-x_eq)**2)

    hist, bins = np.histogram(data, bins = 50)
    p0 = np.sum(hist)

    popt, _ = curve_fit(p_x, bins[1:], hist)
    kappa = popt[0]
    print(kappa)

    hist = [h/p0 for h in hist]
    U = [-kT*np.log(h) for h in hist]

    plt.plot(bins[1:], hist)
    plt.scatter(bins[1:], [p_x(x, kappa)/p0 for x in bins[1:]])
    plt.show()

def PSD(data):
    def lorentz(f, A, B):
        return 1/(A+B*f**2)
    
    def aliased_lorentz(k, dx, c):

        return (dt*dx**2)/(1+c**2-2*c*np.cos(2*np.pi*k))

    freq, Px = periodogram(data, fs = 1/dt)

    def window_inputs(freq, Px, nb=10000):
        steps = int(len(freq)/nb)
        freq = [np.mean(freq[i*steps:(i+1)*steps]) for i in range(nb)]
        Px = [np.mean(Px[i*steps:(i+1)*steps]) for i in range(nb)]
        return freq, Px

    freq, Px = window_inputs(freq, Px)
    N = len(freq)
    indexs = np.array([i/N for i in range(N)])

    plt.plot(freq, Px, c='grey')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e-15, 1e-4])
    

    popt, _ = curve_fit(lorentz, freq, Px)
    A, B = popt[0], popt[1]

    plt.scatter(freq, [lorentz(f, A, B) for f in freq], s = 1, zorder =10)
    freq = list(freq)
    #popt, _ = curve_fit(aliased_lorentz, indexs , Px, bounds = [[0,0], [1,1]])
    #dx, c = popt[0], popt[1]

    #plt.scatter(freq, [aliased_lorentz(i/N, dx, c) for i in range(N)], s = 1)
    
    plt.show()

#potential_well(df['Y-dat'])
PSD(df['X-dat'])

       