import csv
import numpy as np
import matplotlib.pyplot as plt;
from scipy.optimize import curve_fit;
from fitting import fit

#Code for reading and fitting PSD data from FDdat files

#Function reads FDdat file and formats rows into csv cells [Frequency, Px, Py, Sum]
def row_reader(row):
    for num in row:
        lst = (num.split(' '))
        if len(num)==37:
            item1= lst[0]
            item2=lst[1]
        else:
            item3=lst[1]
            item4=lst[2]
    
    new_row =[float(item1), float(item2), float(item3), float(item4)]
    return new_row



#Function reads the name of the file and stores the graph in appropriate folder with useful name.  
def get_label_path(label):
    label = label.split('_')
    if label[0] == 'Sillica':
        material = 'Sillica'
    elif label[0] == 'Vaterite':
        material = 'Vaterite'
    elif label[0] == 'PS':
        material = 'PS'
    else:
        material = 'LC'
    
    power = label[1]
    size = label[2][:-4]
    return material, power, size

#lead path is just the neccesarray jargon to locate directory 
lead_path = '/home/daniel/Documents/dmaciver97/python files/' 
label = 'Vaterite_2.00A_3mu.txt'
material, power, size= get_label_path(label)

path = lead_path+label
        
#convert FDdat file to csv file and then store first 3 coloums as [freq, Px, Py]
with open(path, 'r') as in_file:
    reader = csv.reader(in_file.readlines()[4:])
    
    with open('log.csv', 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerows(row_reader(row) for row in reader)

with open('log.csv', 'r') as file:
    reader = csv.reader(file.readlines())
    freq, Px, Py =[], [], []
    for row in reader:
        freq.append(float(row[0]))
        Px.append(float(row[1]))
        Py.append(float(row[2]))

#Simplified Lorentzian fit of power spectrum with constants A and B
def lorentz(x,A,B):
    return (1/(A+B*x**2))

popt, pcov  = curve_fit(lorentz, freq[5:300], Px[5:300])
A,B = popt
fc = (A/B)**0.5
D = (2*np.pi**2/(B))
Px_model = [D/(2*np.pi**2*(f**2+fc**2))  for f in freq]

fig  = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)  
ax1.plot(freq[1:], Px[1:])
ax1.plot(freq[1:], Px_model[1:], label= f'${np.round(D,3)}/(f^2+{np.round(fc,2)}^2)$')
ax1.legend()
ax1.set_ylabel('Power $[V^2/Hz]$')
ax1.set_xlabel('Frequency [Hz]')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim(1,2.0e03)

popt, pcov  = curve_fit(lorentz, freq[5:300], Py[5:300])
A,B = popt
fc = (A/B)**0.5
D = (2*np.pi**2/(B))
Py_model = [D/(2*np.pi**2*(f**2+fc**2)) for f in freq]

ax2 = fig.add_subplot(1,2,2)
ax2.plot(freq[1:], Py[1:])
ax2.plot(freq[1:], Py_model[1:], label= f'${np.round(D,4)}/(f^2+{np.round(fc,2)}^2)$')
ax2.legend()
#ax2.set_ylabel('Power $[V^2/Hz]$')
ax2.set_xlabel('Frequency [Hz]')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlim(1,2.0e03)
fig.tight_layout()
plt.savefig(f'/home/daniel/Documents/dmaciver97/python files/figures/{material}/{size}@{power}.png')