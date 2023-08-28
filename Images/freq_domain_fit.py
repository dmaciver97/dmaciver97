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
        current_time = datetime()
        self.tag = current_time.now()

    def create_csv_file(self, path):

        def row_reader(self, row):
            line = list(row[0].split())
            new_row = [float(line[0]), float(line[1]), float(line[2]), float(line[3])]

            return new_row
        
        with open(path, 'r') as input:
            reader = csv.reader(input.readlines()[5:])

            with open('tmp_tdat.csv', 'w') as output: 
                writer = csv.writer(output)
                writer.writerows(self.row_reader(row) for row in reader)

        return self.get_data('tmp_fdat.csv')
    
    def get_data(self, file):
        with open(file, 'r') as data:
            reader = csv.reader(data.readlines())
            time, x_QPD, y_QPD = [], [], []
            for row in reader:
                time.append(row[0])
                x_QPD.append(row[1])
                y_QPD.append(row[2])

        return time, x_QPD, y_QPD