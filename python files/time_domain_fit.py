import argparse
import csv
import numpy as np
from scipy.optimize import curve_fit
import scipy
from scipy.signal import correlate
import matplotlib.pyplot as plt
from datetime import datetime

class time_domain:
    ''' 
    Class for processing time domain data from QPD results,
    main two ouputs from this class are the equipartition and 
    auto-correlation fittings. The former providing a rough 
    estimate of the trap stiffness, the latter an estimation
    of the trapping medium's viscosity.

    In order for either method to work we require a conversion
    factor to convert our data from Voltz to m, this conversion 
    factor is computed in the freq_domain.py file. By fitting
    the PSD we can make a comparison between our experimental
    value of D (diffusion constant) and the value of D given
    by the Einstein relation. Both this file and the freq_domain.py
    file should be run in tandem to get accurate results. 
    
    '''
    def __init__(self, beta):
        # Only necessary constant is the conversion factor
        self.beta = beta
        self.kT = 1.38e-23*290
        current_time = datetime()
        self.tag = current_time.now()
    
    def create_csv_file(self, path):
        # 
        with open(path, 'r') as input:
            reader = csv.reader(input.readlines()[5:])

            with open('tmp_tdat.csv', 'w') as output: 
                writer = csv.writer(output)
                writer.writerows(self.row_reader(row) for row in reader)

        return self.get_data('tmp_tdat.csv')
    
    def get_data(self, file):
        with open(file, 'r') as data:
            reader = csv.reader(data.readlines())
            time, x_QPD, y_QPD = [], [], []
            for row in reader:
                time.append(row[0])
                x_QPD.append(row[1]/np.sqrt(self.beta))
                y_QPD.append(row[2]/np.sqrt(self.beta))

        return time, x_QPD, y_QPD

    def row_reader(self, row):
        line = list(row[0].split())
        new_row = [float(line[0]), float(line[1]), float(line[2]), float(line[3])]

        return new_row

    def equipartition_fit(self, data):
        data_eq = np.mean(data)
        variance = [(d-data_eq)**2 for d in data]
        avg_variance = np.mean(variance)
        trap_stiffness = self.kT/avg_variance
        Potential = [0.5*trap_stiffness*v for v in variance]
        plt.close()
        plt.plot([d-data_eq for d in data], Potential,
                 label = r'$\kappa$' + f' = {trap_stiffness} N/m')
        plt.ylabel('Potential [J]')
        plt.xlabel('Displacement [m]')
        plt.legend()
        #plt.savefig(f'python files\Figures\Equipartition\{self.tag}.png')
        plt.show()

        return trap_stiffness
    def get_acf(self, data):
        hatf = np.fft.rfft(data)
        I0 = [f*f.conj() for f in hatf]
        C0 = np.fft.irfft(I0)
        #print(C0[0])
        #C0 = C0 / C0[0]

        return C0

    def get_ACF(self, data, nbatch=50):
        dt = 7.63e-6
        N = len(data)
        n = int(N / nbatch)
        acf_list = []
        for i in range(nbatch):
            tmp = data[i*n:(i+1)*n]
            #print(i*n, (i+1)*n, N)
            avg = np.mean(tmp)
            tmp = [x-avg for x in tmp]
            acf = self.get_acf(tmp)
            acf_list.append(acf)
        t_list = np.arange(0.0, n*dt, dt)
        acf_data = np.array(acf_list)
        acf = [np.mean(acf_data[:,i]) for i in range(n)]
        err = [np.std(acf_data[:,i]) for i in range(n)]
        return acf, err, t_list

    def acf_fit(self, t, A, tau):
        return (A*np.exp(-t/tau))

    def calc_viscosity(self, A, tau):
        return (1.38e-23*290*tau)/(A*6*np.pi*0.785e-6)
    
    def acf_fitting(self, data):
        self.len = 50 
        self.acf_data, err, t_list =  self.get_ACF(data) #acf_data has units m^2, t_list units of seconds 
        self.t_list = t_list.tolist()
        popt, _ = curve_fit(self.acf_fit, self.t_list[:self.len], self.acf_data[:self.len], bounds=(0,1))
        self.A, self.tau = popt #A should have units of m^2, tau units of s 
        model = [self.acf_fit(t, self.A, self.tau) for t in self.t_list]
        mu = self.calc_viscosity(self.A, self.tau)
        self.plot_acf(model, mu)

    #plotting of acf should be done before fitting new data to acf
    def plot_acf(self, model):
        plt.plot(self.t_list[:self.len], self.acf_data[:self.len], label = 'ACF')
        plt.plot(self.t_list[:self.len], model[:self.len], label = f'${np.round(self.A,7)}*exp[-t/{np.round(self.tau,6)}]$')
        plt.legend()
        #plt.savefig(f'python files\Figures\ACF\{self.tag}.png')
        plt.show()
    
if __name__ in '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default=None, help = 'Tdat file to be processed and fitted')
    parser.add_argument('--r', type=float, default=1.57e-6, help='Particle Diamter [m]')
    parser.add_argument('--power', type=float, default=100e-3, help = 'laser power [W]')
    args = parser.parse_args()

