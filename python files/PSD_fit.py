import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt;
from scipy.optimize import curve_fit;

parser = argparse.ArgumentParser()
parser.add_argument('--label', type=str, default='Vaterite_2.00A_3mu.txt', help='FDdat text file')
args = parser.parse_args()

label = args.label

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
    radius  = float(size[:-2])*1e-6
    return material, power, size, radius

#lead path is just the neccesarray jargon to locate directory 
lead_path = '/home/daniel/Documents/dmaciver97/python files/FDdat files/' 
material, power, size, radius= get_label_path(label)

path = lead_path+label
        
#convert FDdat file to csv file and then store first 3 coloums as [freq, Px, Py]
with open(path, 'r') as in_file:
    reader = csv.reader(in_file.readlines()[4:])
    
    with open('log.csv', 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerows(row_reader(row) for row in reader)

with open('log.csv', 'r') as file:
    reader = csv.reader(file.readlines()[:300])
    freq, Px, Py =[], [], []
    for row in reader:
        freq.append(float(row[0]))
        Px.append(float(row[1]))
        Py.append(float(row[2]))

#Simplified Lorentzian fit of power spectrum with constants A and B
def lorentz(x,A,B):
    return (1/(A+B*x**2))

def mod_lorentz(x,a,b,c,fc):
    return (a+(b/((x**2+fc**2)+c**2)))

#Fitting the data to a Lorentzian and then plotting the plots for the X and Y QPD axis
# Calculate conversion factor Beta

gamma0 = 6*np.pi*1.0013e-3*radius 
D_einstein = 1.38e-23*290/gamma0  #units m^2/s

popt, pcov  = curve_fit(lorentz, freq[5:300], Px[5:300])
A,B = popt
fc = (A/B)**0.5
D = (2*np.pi**2/(B)) #units V^2/s
Beta = D/D_einstein #V^2/m^2
Px_model = [D/(Beta*2*np.pi**2*(f**2+fc**2))  for f in freq]
Px = [x/Beta for x in Px]

D_label = str('{:.2e}'.format(D/Beta))
fig  = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)  
ax1.plot(freq[1:], Px[1:])
ax1.plot(freq[1:], Px_model[1:], label= f'${D_label}/(f^2+{np.round(fc,2)}^2)$')
ax1.legend()
ax1.set_ylabel('Power $[m^2/Hz]$')
ax1.set_xlabel('Frequency [Hz]')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim(1,2.0e03)

popt, pcov  = curve_fit(lorentz, freq[5:300], Py[5:300])
A,B = popt
fc = (A/B)**0.5
D = (2*np.pi**2/(B))
Beta = D/D_einstein #V^2/m^2
Py_model = [D/(Beta*2*np.pi**2*(f**2+fc**2)) for f in freq]
Py = [y/Beta for y in Py]

D_label = str('{:.2e}'.format(D/Beta))
ax2 = fig.add_subplot(1,2,2)
ax2.plot(freq[1:], Py[1:])
ax2.plot(freq[1:], Py_model[1:], label= f'${D_label}/(f^2+{np.round(fc,2)}^2)$')
ax2.legend()
#ax2.set_ylabel('Power $[V^2/Hz]$')
ax2.set_xlabel('Frequency [Hz]')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlim(1,2.0e03)
fig.tight_layout()
plt.savefig(f'/home/daniel/Documents/dmaciver97/python files/figures/{material}/{size}@{power}.png')

#popt, pcov =  curve_fit(mod_lorentz, freq[6:300], Px[6:300])
#a,b,c,fc = popt

#model = [a+(b/((x**2+fc**2)+c**2)) for x in Px]

#plt.close('all')
#fig = plt.figure()
#ax = fig.add_subplot()
#ax.plot(freq[1:], Px[1:])
#ax.plot(freq[1:], model[1:])
#ax.set_yscale('log')
#ax.set_xscale('log')
#plt.show()