#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from statsmodels import ccf


def make_ccf(pair1, pair2):
    fwr_diag_cross_product = ccf(pair1, pair2, adjusted=True, fft=True)
    fwr_diag_cross_product /= np.max(abs(fwr_diag_cross_product))
    bck_diag_cross_product = ccf(pair1[::-1], pair2[::-1], adjusted=True, fft=True)[::-1]
    bck_diag_cross_product /= np.max(abs(bck_diag_cross_product))
    return np.r_[bck_diag_cross_product[:-1], fwr_diag_cross_product]

q1_dat, q2_dat, q3_dat, q4_dat = [], [], [], []
t_list = []

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey=True)
plt.suptitle('Quadrant Cross Correlation')
nlags = 50
step_size = 1

tau_ot1 = 3.75e-4
tau_ot2 = 9e-4
with open('./tmp/data_14-May-2024 13_58_16.txt', 'r') as input:
    lines = input.readlines()
    for line in lines[1:-1:step_size]:
        #print(line)
        items =  line.split(',')
        q1_dat.append(float(items[0]))
        q2_dat.append(float(items[1]))
        q3_dat.append(float(items[2]))
        q4_dat.append(float(items[3]))
        t_list.append(float(items[4]))

lag = t_list[2]-t_list[1]
#print('dt = ' + str(lag), 'fs = ' + str(1/lag))

#from statsmodels.graphics.tsaplots import plot_acf, acf

#from test_signal_corr import get_ACF, get_acf

pair1 = [(q1+q2)-(q3+q4) for q1,q2, q3, q4 in zip(q1_dat, q2_dat, q3_dat, q4_dat)]
pair2 = [q1+q2 for q1,q2 in zip(q2_dat, q3_dat)]

ccf_output_diag1 = make_ccf(q1_dat, q4_dat)
ccf_output_diag2 = make_ccf(q2_dat, q3_dat)

pair1 = [q1+q2 for q1,q2 in zip(q1_dat, q2_dat)]
pair2 = [q1+q2 for q1,q2 in zip(q4_dat, q3_dat)]

ccf_output_lr = make_ccf(pair1, pair2)

pair1 = [q1+q2 for q1,q2 in zip(q1_dat, q3_dat)]
pair2 = [q1+q2 for q1,q2 in zip(q2_dat, q4_dat)]

ccf_output_ud = make_ccf(pair1, pair2)

l = len(ccf_output_diag1)//2
if 2*nlags != len(ccf_output_diag1):
    x = np.linspace(0, (nlags+1)*lag, nlags+1)
    ax1.plot(x, ccf_output_diag1[l:l+nlags+1], label ='Left Diagonal')
    ax1.plot(x, ccf_output_diag2[l:l+nlags+1], label ='Right Diagonal')
    ax1.plot(x, ccf_output_lr[l:l+nlags+1], label = 'Left vs Right')
    ax1.plot(x, ccf_output_ud[l:l+nlags+1], label = 'Upper vs Lower')
    #ax1.plot(x, [-np.exp(-lag/tau_ot1) for lag in x], linestyle = '--', c = 'black') 
else:
    x = np.linspace(0, (nlags+1)*lag, nlags)
    ax1.plot(x, ccf_output_diag1[l:l+nlags], label ='Left Diagonal') 
    ax1.plot(x, ccf_output_diag2[l:l+nlags], label ='Right Diagonal') 
    ax1.plot(x, ccf_output_lr[l:l+nlags], label = 'Left vs Right')
    ax1.plot(x, ccf_output_ud[l:l+nlags], label = 'Upper vs Lower')
    #ax1.plot(x, [-np.exp(-lag/tau_ot1) for lag in x], linestyle = '--', c = 'black')  

ax1.legend()
ax1.set_title('Dimer in Circuarly Polarised beam')
ax1.set_xlabel(r'Time Lag $\tau$ [s]')
ax1.set_ylim(-1,1)
ax1.set_ylabel('Quadrant Cross Correlation')

q1_dat, q2_dat, q3_dat, q4_dat = [], [], [], []
t_list = []

with open('./tmp/data_15-May-2024 14_55_24.txt', 'r') as input:
    lines = input.readlines()
    for line in lines[1:-1:step_size]:
        #print(line)
        items =  line.split(',')
        q1_dat.append(float(items[0]))
        q2_dat.append(float(items[1]))
        q3_dat.append(float(items[2]))
        q4_dat.append(float(items[3]))
        t_list.append(float(items[4]))

lag = t_list[2]-t_list[1]

pair1 = [(q1+q2)-(q3+q4) for q1,q2, q3, q4 in zip(q1_dat, q2_dat, q3_dat, q4_dat)]
pair2 = [q1+q2 for q1,q2 in zip(q2_dat, q3_dat)]

ccf_output_diag1 = make_ccf(q1_dat, q4_dat)
ccf_output_diag2 = make_ccf(q2_dat, q3_dat)

pair1 = [q1+q2 for q1,q2 in zip(q1_dat, q2_dat)]
pair2 = [q1+q2 for q1,q2 in zip(q4_dat, q3_dat)]

ccf_output_lr = make_ccf(pair1, pair2)

pair1 = [q1+q2 for q1,q2 in zip(q1_dat, q3_dat)]
pair2 = [q1+q2 for q1,q2 in zip(q2_dat, q4_dat)]

ccf_output_ud = make_ccf(pair1, pair2)

l = len(ccf_output_diag1)//2
if 2*nlags != len(ccf_output_diag1):
    x = np.linspace(0, (nlags+1)*lag, nlags+1)
    ax2.plot(x, ccf_output_diag1[l:l+nlags+1], label ='Left Diagonal')
    ax2.plot(x, ccf_output_diag2[l:l+nlags+1], label ='Right Diagonal')
    ax2.plot(x, ccf_output_lr[l:l+nlags+1], label = 'Left vs Right')
    ax2.plot(x, ccf_output_ud[l:l+nlags+1], label = 'Upper vs Lower')
    #ax2.plot(x, [-np.exp(-lag/tau_ot2) for lag in x], linestyle = '--', c = 'black') 
else:
    x = np.linspace(0, (nlags+1)*lag, nlags)
    ax2.plot(x, ccf_output_diag1[l:l+nlags], label = 'Left Diagonal')
    ax2.plot(x, ccf_output_diag2[l:l+nlags], label = 'Right Diagonal') 
    ax2.plot(x, ccf_output_lr[l:l+nlags], label = 'Left vs Right') 
    ax2.plot(x, ccf_output_ud[l:l+nlags], label = 'Upper vs Lower') 
    #ax2.plot(x, [-np.exp(-lag/tau_ot2) for lag in x], linestyle = '--', c = 'black') 
 
ax2.legend()
ax2.set_ylim(-1, 1)
ax2.set_title('Dimer in Linearly Polarised beam')
ax2.set_xlabel(r'Time Lag $\tau$ [s]')

with open('./tmp/data_15-May-2024 10_44_47.txt', 'r') as input:
    lines = input.readlines()
    for line in lines[1:-1:step_size]:
        #print(line)
        items =  line.split(',')
        q1_dat.append(float(items[0]))
        q2_dat.append(float(items[1]))
        q3_dat.append(float(items[2]))
        q4_dat.append(float(items[3]))
        t_list.append(float(items[4]))

lag = t_list[2]-t_list[1]

pair1 = [q1+q2 for q1,q2 in zip(q1_dat, q4_dat)]
pair2 = [q1+q2 for q1,q2 in zip(q2_dat, q3_dat)]

ccf_output_diag1 = make_ccf(q1_dat, q4_dat)
ccf_output_diag2 = make_ccf(q2_dat, q3_dat)

pair1 = [q1+q2 for q1,q2 in zip(q1_dat, q2_dat)]
pair2 = [q1+q2 for q1,q2 in zip(q4_dat, q3_dat)]

ccf_output_lr = make_ccf(pair1, pair2)

pair1 = [q1+q2 for q1,q2 in zip(q1_dat, q3_dat)]
pair2 = [q1+q2 for q1,q2 in zip(q2_dat, q4_dat)]

ccf_output_ud = make_ccf(pair1, pair2)

l = len(ccf_output_diag1)//2
if 2*nlags != len(ccf_output_diag1):
    x = np.linspace(0, (nlags+1)*lag, nlags+1)
    ax3.plot(x, ccf_output_diag1[l:l+nlags+1], label = 'Left Diagonal')
    ax3.plot(x, ccf_output_diag2[l:l+nlags+1], label = 'Right Diagonal')
    ax3.plot(x, ccf_output_lr[l:l+nlags+1], label = 'Left vs Right')
    ax3.plot(x, ccf_output_ud[l:l+nlags+1], label = 'Upper vs Lower')
else:
    x = np.linspace(0, (nlags+1)*lag, nlags)
    ax3.plot(x, ccf_output_diag1[l:l+nlags], label = 'Left Diagonal') 
    ax3.plot(x, ccf_output_diag2[l:l+nlags], label = 'Right Diagonal') 
    ax3.plot(x, ccf_output_lr[l:l+nlags], label = 'Left vs Right') 
    ax3.plot(x, ccf_output_ud[l:l+nlags], label = 'Upper vs Lower') 

ax3.legend()
ax3.set_ylim(-1, 1)
ax3.set_title('Sphere in Circuarly Polarised beam')
ax3.set_xlabel(r'Time Lag $\tau$ [s]')

plt.show()

