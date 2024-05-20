from PSD_fit import *
from autocorr_fit import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

pathf = '/home/daniel/Documents/dmaciver97/python files/Trapping/FDdat files/'
patht = '/home/daniel/Documents/dmaciver97/python files/Trapping/TDdat files/'

FDfile = open('FDdatfiles.txt', 'r')
TDfile = open('TDdatfiles.txt', 'r')
mu_list = []
acf_list = []
for f, t in zip(FDfile.readlines(), TDfile.readlines()):
    f =  f[:-1]
    t = t[:-1]
    print(f, t)
    psd = PSD_fit(pathf+f)
    auto = autocorr_fit(patht+t)
    freq, x_data, y_data = psd.create_csv_file()
    time, x_sig, s_sig = auto.create_csv_file()
    model, beta = psd.create_fit(freq, x_data)
    psd.plot_fitted_data(freq, x_data, model)
    model_auto, mu = auto.acf_fitting(x_sig, beta)
    acf_list.append(model_auto), mu_list.append(mu)
    auto.plot_acf(model_auto)
FDfile.close()
TDfile.close()
acf =  np.array(acf_list)
acf = [np.mean(acf[:,i]) for i in range(len(acf[0]))]
mu_mean= np.mean(mu_list)

def acf_fit( t, A, tau):
    return (A*np.exp(-t/tau))

def calc_viscosity(A, tau):
        return (1.38e-23*290*tau)/(A*6*np.pi*1.57e-6)
    
t_list = np.arange(0.0, len(acf)*1e-3, 1e-3)

popt, _ = curve_fit(acf_fit, t_list[:2700], acf[:2700])
A, tau = popt
mu = calc_viscosity(A, tau)
print(mu, mu_mean)