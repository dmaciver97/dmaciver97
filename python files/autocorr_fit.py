import argparse
import csv
import numpy as np
from scipy.optimize import curve_fit
import scipy
from scipy.signal import correlate
import matplotlib.pyplot as plt
import statsmodels.api as sm

class autocorr_fit(object):
    def __init__(self, path):
        self.path = path
        self.label = self.path[self.path.rfind("/")+1:]
        #self.material, self.power, self.size, self.radius = self.get_label_path(self.label)

    def row_reader(self, row):
        lst=(row[0].split())
        #print(lst)
        
        new_row = [float(lst[0]), float(lst[1]), float(lst[2]), float(lst[3])]
        #print(new_row)
        return new_row

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
    
    def get_line_count(self, file):
        reader = csv.reader(file.readlines())
        lines = len(list(reader))
        return lines

    def create_csv_file(self):
        with open(self.path, 'r') as in_file:
            reader = csv.reader(in_file.readlines()[5:])
            
            with open('tdat.csv', 'w') as out_file:
                writer = csv.writer(out_file)
                writer.writerows(self.row_reader(row) for row in reader)

        with open('tdat.csv', 'r') as file:
            reader = csv.reader(file.readlines()[:180000])
            #t, x_data, y_data = 0, 0, 0
            t=np.zeros(shape=180000)
            x_data=np.zeros(shape=180000)
            y_data=np.zeros(shape=180000)
            for i, row in enumerate(reader):
                t[i] = row[0]
                x_data[i] = row[1]
                y_data[i] = row[2]

        return t, x_data, y_data
    
    
    def equipartition(self, data):
        data_eq = np.mean(data)
        var = [(d-data_eq)**2 for d in data]
        avg_var = np.mean(var)
        trap_stiffness = 1.38e-23*290/avg_var
        U_x  = [(trap_stiffness/2)*v for v in var]
        print(len(U_x), len(data))
        plt.close()
        plt.plot([d-data_eq for d in data], U_x, label=f'$\kappa_x = {trap_stiffness} N/m$')
        plt.ylabel('Potential [J]')
        plt.xlabel('Displacment [m]')
        plt.legend()
        plt.show()

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
    
    def acf_fitting(self, data, beta):
        self.len = 50 
        data = [x/np.sqrt(beta) for x in data] #Volts to meters
        self.acf_data, err, t_list =  self.get_ACF(data) #acf_data has units m^2, t_list units of seconds 
        self.t_list = t_list.tolist()
        popt, _ = curve_fit(self.acf_fit, self.t_list[:self.len], self.acf_data[:self.len], bounds=(0,1))
        print(popt)
        self.A, self.tau = popt #A should have units of m^2, tau units of s 
        model = [self.acf_fit(t, self.A, self.tau) for t in self.t_list]
        mu = self.calc_viscosity(self.A, self.tau)
        print(mu)
        return model, mu

    #plotting of acf should be done before fitting new data to acf
    def plot_acf(self, model):
        
        plt.plot(self.t_list[:self.len], self.acf_data[:self.len], label = 'ACF')
        plt.plot(self.t_list[:self.len], model[:self.len], label = f'${np.round(self.A,7)}*exp[-t/{np.round(self.tau,6)}]$')
        plt.legend()
        #plt.savefig(f'./figures/Sillica/{self.label}.png')
        plt.show()

if __name__ == "__main__":
    auto = autocorr_fit('/home/daniel/Documents/dmaciver97/python files/Trapping/TDdat files/run9.TDdat')
    time, x, y = auto.create_csv_file()
    #data = [np.sqrt(i**2+j**2) for i, j in zip(x, y)]
    model, mu = auto.acf_fitting(x, 1e5)
    auto.plot_acf(model)