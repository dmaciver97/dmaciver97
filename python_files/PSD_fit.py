import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt;
from scipy.optimize import curve_fit
from fitnest import *

class PSD_fit(object):
    def __init__(self, path):
        self.path = path
        self.label = self.path[self.path.rfind("/")+1:]
        #self.material, self.power, self.size, self.radius = self.get_label_path(self.label)

#Function reads FDdat file and formats rows into csv cells [Frequency, Px, Py, Sum]
    def row_reader(self, row):
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
    def get_label_path(self, label):
        label = label.split('_')
        if label[0] == 'Sillica':
            self.material = 'Sillica'
        elif label[0] == 'Vaterite':
            self.material = 'Vaterite'
        elif label[0] == 'PS':
            self.material = 'PS'
        else:
            self.material = 'LC'
        
        self.power = label[1]
        self. size = label[2][:-4]
        self.radius  = float(self.size[:-2])*1e-6
        return self.material, self.power, self.size, self.radius

#convert FDdat file to csv file and then store first 3 coloums as [freq, Px, Py]
    def create_csv_file(self):
        with open(self.path, 'r') as in_file:
            reader = csv.reader(in_file.readlines()[4:])
            
            with open('log.csv', 'w') as out_file:
                writer = csv.writer(out_file)
                writer.writerows(self.row_reader(row) for row in reader)

        with open('log.csv', 'r') as file:
            reader = csv.reader(file.readlines()[:300])
            col1, col2, col3 =[], [], []
            for row in reader:
                print(row)
                col1.append(float(row[0]))
                col2.append(float(row[1]))
                col3.append(float(row[2]))

            self.fitter = fit(['A', 'B'], col1[4:280], col2[4:280])
        return col1, col2, col3

#Simplified Lorentzian fit of power spectrum with constants A and B
    def lorentz(self,x,A,B):
        return (1/(A+B*x**2))
#Fitting the data to a Lorentzian and then plotting the plots for the X and Y QPD axis
# Calculate conversion factor Beta

    def create_fit(self, col1, col2):
        gamma0 = 6*np.pi*1.0013e-3*1.57e-6
        D_einstein = 1.38e-23*290/gamma0  #units m^2/s
        self.fitter.run()

        file = open('./results/info/post_summary.csv')

        reader = csv.reader(file)
        headings = next(reader)
        row1 = next(reader)
        A=float(row1[0])
        B = float(row1[5])
        
        file.close()

        self.fc = (A/B)**0.5
        self.D = (2*np.pi**2/(0.9*B)) #units V^2/s
        Beta = self.D/D_einstein #V^2/m^2
        print(np.sqrt(Beta))
        model = [self.D/(np.pi**2*(f**2+self.fc**2))  for f in col1]
        return model, Beta

    def plot_fitted_data(self, col1, col2, model):
        D_label = str('{:.2e}'.format(self.D))
        fig  = plt.figure()
        ax= fig.add_subplot() 
        ax.plot(col1[:280], col2[:280])
        ax.plot(col1[:280], model[:280], label= f'${D_label}/(f^2+{np.round(self.fc,2)}^2)$')
        ax.legend()
        ax.set_ylabel('Power $[m^2/Hz]$')
        ax.set_xlabel('Frequency [Hz]')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim(10,2.0e03)
        plt.savefig(f'./figures/Sillica/{self.label}.png')
        plt.close()