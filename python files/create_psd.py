#!/usr/bin/python3

from ast import alias
from matplotlib.mlab import psd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy.signal import periodogram
from scipy.optimize import curve_fit, least_squares

def create_psd(data_file):
        
    x_dat = []
    y_dat = []
    z_dat = []
    t_list = []

    if data_file[-3:] == 'txt':
        print('Loading QPD data')
        test = data_file.split('/')
        tag = test[2][:-3]
        with open(data_file, 'r') as input:
            lines = input.readlines()
            for line in lines[1:-1]:
                items =  line.split(',')
                x_dat.append(float(items[0]))
                y_dat.append(float(items[1]))
                z_dat.append(float(items[2]))
                t_list.append(float(items[3]))
        
        dt = t_list[10]-t_list[9]
        fs = 1/dt

    elif data_file[-3:] == 'npy':
        print('Loading trajectory file')
        test = data_file.split('/')
        tag = test[2][:-8]
        traj = np.load(data_file)
        dt = 1e-5
        for i, frame in enumerate(traj[1000:]):
            x_dat.append(frame[0])
            y_dat.append(frame[1])
            t_list.append(dt*i)
        
        fs = 1/dt

    def loretzian (x, A, B):
            return 1/(A+B*x**2)
    
    def alias_loretzian(k, dx, c, i):
            N = len(freq_block[i:])
            return (dt*dx**2)/(1+c**2-2*c*np.cos(np.pi*(k+1)/N))
    
    def res_lorentz(x, k, y, i):
        y_ex = alias_loretzian(k, x[0], x[1], i)
        sig = y_ex**-1/(nb)**0.5
        return ((y**-1 - y_ex**-1)/sig)

    def itereative_fitting(freq, psd):
        '''
        Fit Loretzian while shifting nb iteratively 
        in order to find the range that minimises the 
        square residuals, then plot the result. 
        '''
       
        x0 = np.array([1,1])
        errors = []

        for i in range(5):
            continue
            N = len(freq[i:])
            K = np.arange(0, N, 1)
            result = least_squares(res_lorentz, x0, args=(K, np.array(psd[i:]), i))
            errors.append(result.cost)

        index = 1#errors.index(min(errors))

        N = len(freq[index:])
        K = np.arange(0, N, 1)
        result = least_squares(res_lorentz, x0, args=(K, np.array(psd[index:]), index))

        return index, result.x[0], result.x[1]
    
    fig, ax1 = plt.subplots()
    #ax2 = fig.add_axes([0.22, 0.2, 0.4, 0.3])

    nb = 500

    N = len(x_dat)

    xdft = fft(x_dat)
    xdft = xdft[:int(N/2)+1]
    n = len(xdft)
    psd_x = [1/(fs*N)*abs(x)**2 for x in xdft]
    psd_block_x = [np.mean(psd_x[i*nb:(i+1)*nb]) for i in range(0, int(n/nb))]

    freq = np.linspace(0, int(fs/2), int(len(psd_x)))
    freq_block = [np.mean(freq[i*nb:(i+1)*nb]) for i in range(0, int(n/nb))]
    freq_block = [x for x in freq_block if str(x) != 'nan']
    psd_block_x = psd_block_x[:len(freq_block)]
    
    #ax1.plot(freq, psd_x)
    index, dx_x, c_x = itereative_fitting(freq_block, psd_block_x)
    #popt, _ = curve_fit(alias_loretzian, np.arange(0,len(freq_block[index:]), 1), 
    #                    psd_block_x[index:])
    #corner_freq_x = (dx_x/c_x)**0.5
    #D_x = 2*np.pi**2/c_x
    #dx_x, c_x = popt[0], popt[1]

    corner_freq_x = -np.log(c_x)*fs/(2*np.pi)
    D_x = (dx_x**2)*2*np.pi*corner_freq_x/(1-c_x**2)

    ax1.scatter(freq_block[index:], psd_block_x[index:], 
                zorder=10, s = 3, label = r'$P_{x,blocked}$')
    ax1.plot(freq_block[index:], [alias_loretzian(k, dx_x, c_x, index) for k in np.arange(0,len(freq_block[index:]), 1)], 
             zorder=5, label = r'$f_{c,x}$' + f'={abs(corner_freq_x):.3g}')
    #ax1.plot(freq_block[index:], [loretzian(f, dx_x, c_x) for f in freq_block[index:]], 
    #         zorder=10, label = r'$f_{c,x}$' + f'={abs(corner_freq_x):.3g}')


    ydft = fft(y_dat)
    ydft = ydft[:int(N/2)+1]
    n = len(ydft)
    psd_y = [1/(fs*N)*abs(y)**2 for y in ydft]
    psd_block_y = [np.mean(psd_y[i*nb:(i+1)*nb]) for i in range(0, int(n/nb))]

    #ax1.plot(freq, psd_y)
    index, dx_y, c_y = itereative_fitting(freq_block, psd_block_y)
    #popt, _ = curve_fit(alias_loretzian, np.arange(0,len(freq_block[index:]), 1),
    #                    psd_block_y[index:])
    #corner_freq_y = (dx_y/c_y)**0.5
    #D_y = 2*np.pi**2/c_y
    #dx_y, c_y = popt[0], popt[1]

    corner_freq_y = -np.log(c_y)*fs/(2*np.pi)
    D_y = (dx_y**2)*2*np.pi*corner_freq_x/(1-c_y**2)

    ax1.scatter(freq_block[index:], psd_block_y[index:], 
                zorder=10, s = 3 ,label = r'$P_{y, blocked}$')
    ax1.plot(freq_block[index:], [alias_loretzian(k, dx_y, c_y, index) for k in np.arange(0,len(freq_block[index:]), 1)], 
             zorder=5, label = r'$f_{c,y}$' + f'={abs(corner_freq_y):.3g}')
    #ax1.plot(freq_block[index:], [loretzian(f, dx_y, c_y) for f in freq_block[index:]], 
    #         zorder=10, label = r'$f_{c,x}$' + f'={abs(corner_freq_y):.3g}')

    

    #psd_xy = [np.real((p_x**0.5)*np.conj(p_y)**0.5) for p_x, p_y in zip(psd_block_x, psd_block_y)]
    #pxy_px_py = [p_xy/(p_x*p_y)**0.5 for p_xy, p_x, p_y in zip(psd_xy, psd_block_x, psd_block_y) ]
    

    gamma = 6*np.pi*1e-3*1e-6
    
    print(abs(corner_freq_x), abs(corner_freq_x)*2*np.pi*(1.38e-23*298/D_x))
    print(abs(corner_freq_y), abs(corner_freq_y)*2*np.pi*(1.38e-23*298/D_y))
    
    #ax1.plot(freq_block[1:], [alias_loretzian(k, dx_y, c_y) for k in K[:]], zorder=10, label = r'$f_{c,y}$' + f'={abs(corner_freq_y):.3g}')
    #ax1.plot(freq_block, [D_y/(2*np.pi**2*(corner_freq_y**2+f**2)) for f in freq_block], zorder=10)

    #ax1.scatter(corner_freq_x, alias_loretzian(abs(corner_freq_x), dx_x, c_x), c='g', zorder=15, s = 50, label = r'$f_{c,x}$' + f'={abs(corner_freq_x):.3g}')
    #ax1.scatter(corner_freq_y, alias_loretzian(abs(corner_freq_y), dx_y, c_y), c='r', zorder=15, s = 50, label = r'$f_{c,y}$' + f'={abs(corner_freq_y):.3g}')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.legend(loc = 'best')
    ax1.set_xlabel('Frequency [Hz]')
    if data_file[-3:] == 'txt':
        ax1.set_ylabel('P(f) $[V^2/Hz]$')
    elif data_file[-3:] == 'npy':
        ax1.set_ylabel('P(f) $[m^2/Hz]$')
    
    #ax1.set_xlim(10, 1e5)

    #ax2.scatter(freq_block, pxy_px_py)
    #ax2.set_ylabel(r'$P_{x,y}/P_xP_y$')
    #ax2.set_yscale('log')
    plt.savefig(f'./output/PSDs/psd_from_[{tag}].png')
    plt.close('all')

    return abs(corner_freq_x), abs(corner_freq_y)

if __name__ in '__main__':
    create_psd('QPD_data_21-02-24.txt')