import argparse
import csv
import numpy as np
from scipy.optimize import curve_fit
import scipy
from scipy.signal import correlate
import matplotlib.pyplot as plt
from datetime import datetime

class freq_domain:
    '''
    Class for running frequency domain analysis from QPD data, 
    '''
    def __init__(self):
        self.kT = 1.38e-23*290
        current_time = datetime.now()
        self.tag = current_time.now()

    def create_csv_file(self, file):

        def row_reader(row):
            line = list(row[0].split())
            line += list(row[1].split())
            new_row = [float(line[0]), float(line[1]), float(line[2]), float(line[3])]

            return new_row
        main_path = 'C:\\Users\\xbb16146\\dmaciver97\\data\\Trapping\\'
        path = main_path+file
        with open(path, 'r') as input:
            reader = csv.reader(input.readlines()[5:])

            with open('tmp_fdat.csv', 'w', newline='') as output: 
                writer = csv.writer(output)
                writer.writerows(row_reader(row) for row in reader)

        return self.get_data('.\\tmp_fdat.csv')
    
    def get_data(self, file):
        with open(file, 'r') as data:
            reader = csv.reader(data.readlines())
            freq, x_QPD, y_QPD = [], [], []
            for row in reader:
                freq.append(float(row[0]))
                x_QPD.append(float(row[1]))
                y_QPD.append(float(row[2]))

        return freq, x_QPD, y_QPD
    
    def plot_psds(self, freq, x_fdat, y_fdat):
        plt.plot(freq, x_fdat)
        plt.plot(freq, y_fdat)
        plt.xscale('log')
        plt.yscale('log')
        plt.show()

    def fit_psd(self, freq, fdat):
        
        def simplified_psd(f, A, B):
            return 1/(A+B*f**2)
        
        def evaluate_fit(fit, fdat):
            SS_res = np.sum([(y-f)**2 for y,f in zip(fdat, fit)])
            SS_tot = np.sum([(y-np.mean(fdat))**2 for y in fdat])
            return 1-(SS_res/SS_tot)
        
        r2 = 0
        for i in range(5):
            popt, _ = curve_fit(simplified_psd, freq[i:], fdat[i:])
            A, B = popt
            fit  = [simplified_psd(f, A, B) for f in freq[i:]]
            r2_new = evaluate_fit(fit, fdat[i:])
            print(r2_new)
            if r2_new > r2:
                r2 = r2_new
                tmp = i
        print(tmp)
        plt.plot(freq[tmp:], [simplified_psd(f, A, B) for f in freq[tmp:]])
        plt.plot(freq[tmp:], fdat[tmp:])
        plt.xscale('log')
        plt.yscale('log')
        plt.show()


if __name__ in '__main__':
    sim = freq_domain()
    freq, x_QPD, y_QPD = sim.create_csv_file('2021\\02-12-2021\\Sillica+Water_600W_157.FDdat')
    sim.fit_psd(freq, x_QPD)