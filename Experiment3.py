#!/usr/bin/python3

import heapq
import math
import csv
import numpy as np
import pickle
import pylab as plt
from statistics import mean, stdev
from scipy.interpolate import splrep, splev, splint
from scipy.integrate import tplquad
from itertools import product
from numpy.linalg import norm
import random
import argparse
import itertools
import xml.etree.cElementTree as ET
from xml.dom import minidom
import multiprocessing


from Beam import *
from MSTM import *
from Particle import *

neighbor_list = [np.array(vec) for vec in product([-1, 0, 1], repeat=3)]




################################################################################
################################################################################
################################################################################
class Experiment3(object):


    def __init__(self, pklfile='fit.pkl', angle_list=None):
        
        self.theta_list = None
        self.n_list = self.set_measurement_angles        
        self.sigma = 0.1e1
        self.fit_dict = self.read_pickle(pklfile)
        self.N = 10000
        self.mstm = MSTM()

        if (angle_list is None):
            angle_list = [15.0, 45.0, 90.0]
        self.set_measurement_angles(angle_list)

        
################################################################################
    def read_pickle(self, filename):
        pklfile = open(filename, 'rb')
        fit_dict = pickle.load(pklfile)
        pklfile.close()
        return fit_dict


################################################################################
    def set_measurement_angles(self, theta_list):

        self.theta_list = theta_list
        
        n_list = []
        Iraw_measurements = [[] for theta in theta_list]
        Deriv_measurements =[[] for theta in theta_list]
        sum_prob = 0.0
        for n, fit in self.fit_dict.items():
            Iraw = [ float(splev(theta, fit) ) for theta in theta_list]
            Deriv = [float(splev(theta, fit, der =1)) for theta in theta_list]
            for index, I in enumerate(Iraw):
                Iraw_measurements[index].append(I)
            prob = 1.0
            sum_prob += prob
            tmp_dict = {'n': n, 'I_raw': Iraw, 'sigma': self.sigma, 'p': sum_prob}
            n_list.append(tmp_dict)
        
        I_mean = [mean(Iraw)  for Iraw in Iraw_measurements]
        I_std  = [stdev(Iraw) for Iraw in Iraw_measurements]
        Deriv_mean = [mean(Iraw)  for Iraw in Iraw_measurements]
        Deriv_std  = [stdev(Iraw) for Iraw in Iraw_measurements]
        self.I_avg = I_mean
        self.I_std = I_std
        
        I_dat = []
        for Iraw_list, mu, sig in zip(Iraw_measurements, I_mean, I_std):
            I_list = [(Iraw-mu)/sig for Iraw in Iraw_list]
            I_dat.append(I_list)
#            print(f'  mu={mean(I_list)}  std={stdev(I_list)}')
        self.I_max  = [max(I_list)   for I_list in I_dat]
        self.I_min  = [min(I_list)   for I_list in I_dat]

        
        #print(I_dat)
        self.y_list = []
        for n in n_list:
            n['y'] = [(Iraw-mu)/sig for Iraw, mu, sig in zip(n['I_raw'], I_mean, I_std)]
            n['p'] = n['p'] / sum_prob
            self.y_list.append(n['y'])
        self.n_list = n_list.copy()

        self.plot_nvsy()

        print(self.dy)
        
        return n_list
###############################################################################
###############################################################################
    def plot_nvsy(self):
        
        dist_list = []
        for I in self.y_list:
            for  y in self.y_list:
                dist = np.sqrt((y[0]-I[0])**2+(y[1]-I[1])**2+(y[2]-I[2])**2)
                if dist !=0:
                    dist_list.append(dist)
        self.dy =min(dist_list)
    
        dot_list = []
        n1=[0.00000, 0.00000, 1.0000]
        for dot_dict in self.n_list:    
            n2 = dot_dict['n']
            dot  = np.dot(n1,n2)
            if dot != 0:
                dot_list.append(dot)


        #fig = plt.figure()
        #ax = fig.add_subplot()
        #ax.scatter(dist_list, dot_list)
        #for label, dist in enumerate(dist_list):
            #ax.text(dist_list[label], dot_list[label], str(label+1))

        #ax.set_xlabel('Difference in signals')
        #ax.set_ylabel('Difference in orientation')
        #ax.set_title('$d\hat{n}$ vs $dy$ for $\hat{n}_1$')
        
        #plt.savefig('/home/daniel/programs/pyMSTM/tests/twoparticle/dir4/PhD_Notes/Plots/RCE files/dnvsdy.png')  
        #plt.show()

################################################################################
    def get_y(self, I_raw):

        I_scaled = [(raw-avg)/std for raw, avg, std in zip(I_raw, self.I_avg, self.I_std)]
        y = np.array(I_scaled)

        return y
        

################################################################################
    def get_dist(self, s1, s2):
        dot = np.dot(np.array(s1), np.array(s2))
        #print(dot)
        dist = np.arccos(dot) 
        return dot

        
################################################################################
    def get_distances(self, s):

        s_vec = np.array(s)
        dist_list = []
        
        for n in self.n_list:
           
            n_vec = np.array(n['n'])
           
            dot = np.dot(s_vec, n_vec)
            dist = np.arccos(dot)
            dist_list.append(dist)
#            print(n['n'], dist)
        dist_min = min(dist_list)
        index = dist_list.index(dist_min)
#        print(index, dist_min)
        return index

###############################################################################
################################################################################
    def estimate_orientation(self, I_raw, s, epsilon):
    #    print('---')
        prob_list = []

        sigma1, sigma2, sigma3 = [epsilon*I for I in self.I_avg]

        sigma1 /=self.I_std[0]
        sigma2 /=self.I_std[1]
        sigma3 /=self.I_std[2]
        #print(sigma1, sigma2, sigma3)
        exp_list = []

        for n_dict in self.n_list:

            n = n_dict['n']
            dot = np.dot(s,n)
            exp_list.append(np.exp((dot-1)))

        dot = []
        
        for count, n_dict in enumerate(self.n_list):

            n = n_dict['n']
            #print(n)
            boltz = exp_list[count]/np.sum(exp_list)
            dot.append(boltz)
            #print(dot)
            #dot.append(np.dot(s,n))
            #print(dot_norm)
            I_scaled = [(raw-avg)/std for raw, avg, std in zip(I_raw, self.I_avg, self.I_std)]
            I_scaled = np.array(I_scaled)
           
            y_ref = np.array( n_dict['y'] )
            #print(I_scaled, y_ref)
            I1 = I_scaled[0] + random.uniform(-1,1)*sigma1
            I2 = I_scaled[1] + random.uniform(-1,1)*sigma2
            I3 = I_scaled[2] + random.uniform(-1,1)*sigma3
           
            gauss1 = np.exp(-(I1-y_ref[0])**2/2*sigma1**2)/(np.sqrt(2*np.pi*sigma1**2))**3
            gauss2 = np.exp(-(I2-y_ref[1])**2/2*sigma2**2)/(np.sqrt(2*np.pi*sigma2**2))**3
            gauss3 = np.exp(-(I3-y_ref[2])**2/2*sigma3**2)/(np.sqrt(2*np.pi*sigma3**2))**3
            #print(gauss1, gauss2, gauss3)
            probyn = (gauss1*gauss2*gauss3)
            #print(probyn)
            prob = probyn*boltz
            prob_list.append(prob)
            #print(prob_list)
            #print(n, prob)
            #print(f'  {s}: {norm(I_scaled-I_ref)}')

        plt.hist(dot,bins = 10,  density =True)
        plt.show()
        
        py = np.sum(prob_list)
        for x in range(len(prob_list)):
            prob_list[x] = prob_list[x]/py
        prob_max = max(prob_list)
        #print(np.sum(prob_list))
        #print('-----')
        #print(proby)

        #print(self.score_run(dot, prob_list))
      
        index = prob_list.index(prob_max)
        #print(index, prob_max)
        return index, prob_list


##############################################################################
    def write_xml(self, xml, outfile):

        raw_str = ET.tostring(xml)
        refined_xml = minidom.parseString(raw_str)
    #    print(refined_xml.toprettyxml(indent="  "))

        f = open(outfile, 'w')
        f.write(refined_xml.toprettyxml(indent="  "))
        f.close()

###############################################################################
    def score_run(self, dot, estimate):
        print(dot,estimate)
        K_l = 0
        angle_list = self.theta_list
        
        for x in range(len(dot)):
            if dot[x] <= 0:
                div = - 10

                div = -estimate[x]*np.log(dot[x]/estimate[x])
            print(div)
            epsilon = abs(dot[x]-estimate[x])
            K_l += div
    

        return K_l


################################################################################
    def get_slist(self):
        '''This function generates a set of reference orientations.'''

        s0_list = []

        s0_list.append( np.array([0.2958759, 0.2958759, 0.9082483]) )
        s0_list.append( np.array([0.9082483, 0.2958759, 0.2958759]) )
        s0_list.append( np.array([0.2958759, 0.9082483, 0.2958759]) )
        s0_list.append( np.array([1.0, 0.0, 0.0]) )
        s0_list.append( np.array([0.0, 1.0, 0.0]) )
        s0_list.append( np.array([0.0, 0.0, 1.0]) )

        s_list = []
        for f_vec in itertools.product([1.0, -1.0], repeat=3):
            #print(f_vec)
            for s0_vec in s0_list:
                s_vec = np.array( [f*s0 for f, s0 in zip(f_vec, s0_vec)] )
                s_list.append(s_vec)

        return s_list

################################################################################
    def get_scatter_fit(self, s):
        pos1 = list(  s * radius1/l)
        pos2 = list(- s * radius2/l)
        cluster.particle_list[0].set_position(pos1)
        cluster.particle_list[1].set_position(pos2)

        cluster.exec()
        theta_dat, phi_dat, I_dat = cluster.get_scatter(normalized=False, degrees=True)
        fit = splrep(theta_dat, I_dat)

        return fit

################################################################################
    def get_data(self, s_list):

        xml = ET.Element('Simulation')

        fit_dict = {}
        for index, s in enumerate(s_list):
    #        print(s)
            s_xml = ET.SubElement(xml, 's', {f's{i}': f'{x}' for i,x in enumerate(s)})
            pos1 = list(  s * radius1/l)
            pos2 = list(- s * radius2/l)
            cluster.particle_list[0].set_position(pos1)
            cluster.particle_list[1].set_position(pos2)

            cluster.exec()
            theta_dat, phi_dat, I_dat = cluster.get_scatter(normalized=False, degrees=True)

            for theta, phi, I in zip(theta_dat, phi_dat, I_dat):
                tmp_dict = {'theta': f'{theta}', 'phi': f'{phi}',  'I': f'{I[0]}'}
                data_xml = ET.SubElement(s_xml, 'data', tmp_dict)
            fit = splrep(theta_dat, I_dat)
            print(fit)
            fit_dict[tuple(s)] = fit
            plt.plot(theta_dat, I_dat)

        ymin = 1.0e-5
        ymax = 1.0e2
        
        for angle in angle_list:
            plt.vlines(angle, ymin, ymax, color='black', ls='dotted')
        plt.vlines(45.0, ymin, ymax, color='black', ls='dotted')
        plt.vlines(90.0, ymin, ymax, color='black', ls='dotted')
        plt.xlabel(r'$\theta$')
        plt.ylabel(r'$I(\theta)$')
        plt.yscale('log')
        plt.xlim([0.0, 180.0])
        plt.ylim([ymin, ymax])
        plt.savefig('scatter_reference.png')
        plt.close('all')
        plt.show()

        return fit_dict, xml

###############################################################################
################################################################################
    def set_integration_list(self):

        self.box = 5.0**self.sigma
        cell_list = []
        for s in self.s_list:
            cell = [math.floor(x/self.box) for x in s['s']]
            cell_list.append(cell)
        
        integration_list = []
        for cell in cell_list:
            tmp = np.array(cell, dtype=int)
            for neighbor in neighbor_list:
                integration_list.append(tuple(tmp+neighbor))
        self.integration_list = set(integration_list)
#        for cell in integration_list:
#            print(cell)
#            print(f' {cell[0]*self.sigma}, {(cell[0]+1)*self.sigma}')

        
    
################################################################################
    def I_response(self, angle1, angle2, angle3, fit_dict):
        angle_list = [angle1, angle2, angle3]

        I1_raw = []
        I2_raw = []
        I3_raw = []
        for s, fit in fit_dict.items():
            #        I1 = float( splev(angle1, fit) )
            #        print(euler_angles, type(I1), I1)
            I1_raw.append(float(splev(angle1, fit)))
            I2_raw.append(float(splev(angle2, fit)))
            I3_raw.append(float(splev(angle3, fit)))
    
        I1_avg, I1_std = mean(I1_raw), stdev(I1_raw)
        I2_avg, I2_std = mean(I2_raw), stdev(I2_raw)
        I3_avg, I3_std = mean(I3_raw), stdev(I3_raw)
        I1_dat = [(I-I1_avg)/I1_std for I in I1_raw]
        I2_dat = [(I-I2_avg)/I2_std for I in I2_raw]
        I3_dat = [(I-I3_avg)/I3_std for I in I3_raw]

        xml = ET.Element('Response')
        probe1_xml = ET.SubElement(xml, 'Probe',
                                   {'id': f'{0}',
                                    'theta': f'{angle1}',
                                    'I_avg': f'{I1_avg}',
                                    'I_std': f'{I1_std}'})
        probe2_xml = ET.SubElement(xml, 'Probe',
                                   {'id': f'{1}',
                                    'theta': f'{angle2}',
                                    'I_avg': f'{I2_avg}',
                                    'I_std': f'{I2_std}'})
        probe3_xml = ET.SubElement(xml, 'Probe',
                                   {'id': f'{2}',
                                    'theta': f'{angle3}',
                                    'I_avg': f'{I3_avg}',
                                    'I_std': f'{I3_std}'})
        data_xml = ET.SubElement(xml, 'Data')

        for s, Iraw1, Iraw2, Iraw3, I1, I2, I3 in zip(fit_dict.keys(), I1_raw, I2_raw, I3_raw, I1_dat, I2_dat, I3_dat):
            tmp = {f's{i}': f'{x}' for i,x in enumerate(s)}
            s_xml = ET.SubElement(data_xml, 's', tmp)
            raw_xml = ET.SubElement(s_xml, 'raw',
                                       {'I0': f'{Iraw1}', 'I1': f'{Iraw2}', 'I2': f'{Iraw3}'})
            scaled_xml = ET.SubElement(s_xml, 'scaled',
                                       {'I0': f'{I1}', 'I1': f'{I2}', 'I2': f'{I3}'})

        write_xml(xml, 'fit2.xml')    



        return I1_dat, I2_dat, I3_dat

################################################################################
    def pq(self, I):
    
        f = 0.0
        for s in self.n_list:
            dI = I - s['y']
            dI2 = np.dot(dI, dI)
            sig2 = s['sigma']**2
            f += np.exp(-0.5*dI2/sig2)
        f /= np.sqrt(2.0*np.pi*sig2)**3
        
        return f


################################################################################
    def Hq_int(self, I0, I1, I2):
    
        I = np.array([I0, I1, I2])
        p = self.pq(I)
        H = 0.0
        if (p > 0.0):
            H = p * np.log(p)
            
        return H

    
################################################################################
    def get_Hq(self, theta_list):
        
        self.set_measurement_angles(theta_list)
        Hq = tplquad(self.Hq_int,
                     self.I_min[0], self.I_max[0],
                     self.I_min[1], self.I_max[1],
                     self.I_min[2], self.I_max[2])
        print(theta_list, Hq)
        return Hq[0]


################################################################################
    def ln_pq(self, I):

        sig2 = self.sigma**2
        f = 0.0
        for s in self.s_list:
            dI = I - s['I']
            dI2 = np.dot(dI, dI)
            f += np.exp(-0.5*dI2/sig2)
        
        return np.log(f) - 1.5*np.log(2.0*np.pi*sig2)


################################################################################
    def get_Hq_MC(self, theta_list):

        self.set_measurement_angles(theta_list)
        count = 0
        Hq = 0.0
#        for i in range(self.n*len(self.s_list)):
#            prob = random.random()
#            for s_tmp in self.s_list:
#                if (prob < s_tmp['p']):
#                    s = s_tmp
#                    break
        for s in self.s_list:
            dI_list = np.random.normal(0.0, self.sigma, (self.N, 3))
            for dI in dI_list:
                I = s['I'] + dI
                Hq += self.ln_pq(I)
                count += 1
        Hq /= float(count)
        print(theta_list, Hq)
        return Hq


################################################################################
    def get_Hq_cell(self, theta_list):

        self.set_measurement_angles(theta_list)
        print(len(self.integration_list))
        Hq = 0.0
        for cell in self.integration_list:
            print(cell)
            tmp = tplquad(self.Hq_int,
                         cell[0]*self.box, (cell[0]+1)*self.box,
                         cell[1]*self.box, (cell[1]+1)*self.box,
                         cell[2]*self.box, (cell[2]+1)*self.box)
            Hq += tmp[0]
        print(theta_list, Hq)
        return Hq

    
################################################################################
################################################################################
################################################################################
if __name__ == '__main__':

    from subprocess import check_call
    import quaternion
    from quaternion import as_rotation_matrix
    
    a = 1.0e3   # sphere I radius / nm
    lam = 2.0
    logfile = open('junk.log', 'w')
    
    # determine scattering patterns for multiple orientations
    cmd  = ['./test_orientation2.py']
    cmd += ['--pkl', f'fit.5.pkl']
    cmd += ['--R1', f'{a}']
    cmd += ['--R2', f'{a/lam}']
    cmd += ['--z1', f'{a}']
    cmd += ['--z2', f'{-a/lam}']
    print(cmd)
    check_call(cmd, stdout=logfile, stderr=logfile)

    
    sim = Experiment3('fit.5.pkl')
#    print(sim.fit_dict.keys())

    theta1, theta2, theta3 = 10.0, 45.0, 90.0
    theta_list = [theta1, theta2, theta3]
    sigma = 0.1e1
    n_list = sim.set_measurement_angles(theta_list)    

#    for n in n_list:
#        print(n)

    traj_file = 'traj.csv.npy'
    traj = np.load(traj_file)
    frame = traj[1000]
    r = np.array([frame[0], frame[1], frame[2]])
    q = np.quaternion(frame[3], frame[4], frame[5], frame[6])
    R = as_rotation_matrix(q)
    s = np.array([R[0,2], R[1,2], R[2,2]])
    print(s)
    
#    for count, frame in enumerate(traj):
#        r = np.array([frame[0], frame[1], frame[2]])
#        q = np.quaternion(frame[3], frame[4], frame[5], frame[6])
#        R = as_rotation_matrix(q)
#        s = np.array([R[0,2], R[1,2], R[2,2]])



#    Hq = sim.Hq_int(0.1, 0.1, -0.5)
#    Hq = sim.get_Hq(theta_list)
#    print(Hq)
    
    
################################################################################
################################################################################
################################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--lam', type=float, default=1064.0, help='(vacuum) wavelength of incident beam / nm')
parser.add_argument('--CB', type=float, default=0.0, help='inverse beam width')
parser.add_argument('--R1', type=float, default=1000.0, help='radius of particle 1 / nm')
parser.add_argument('--R2', type=float, default=500.0, help='radius of particle 2 / nm')
parser.add_argument('--n1', type=complex, default=1.59+0.0j, help='refractive index of particle 1')
parser.add_argument('--n2', type=complex, default=1.59+0.0j, help='refractive index of particle 2')
parser.add_argument('--n', type=complex, default=1.33+0.0j, help='refractive index of medium')
parser.add_argument('--x1', type=float, default=0.0, help='x-position of sphere 1 / nm')
parser.add_argument('--y1', type=float, default=0.0, help='y-position of sphere 1 / nm')
parser.add_argument('--z1', type=float, default=-1000.0, help='z-position of sphere 1 / nm')
parser.add_argument('--x2', type=float, default=0.0, help='x-position of sphere 2 / nm')
parser.add_argument('--y2', type=float, default=0.0, help='y-position of sphere 2 / nm')
parser.add_argument('--z2', type=float, default=500.0, help='z-position of sphere 2 / nm')
parser.add_argument('--stokes', type=float, nargs='+', default=[1.0, 0.0, 0.0, 0.0], help='Stokes parameter of incoming beam')
parser.add_argument('--angle_list', type=float, default=[10.0, 45.0, 90.0])

args = parser.parse_args()

CB = args.CB

angle_list = args.angle_list

stokes = args.stokes

n1, n2, n = args.n1, args.n2, args.n 

lam = args.lam
k = 2*np.pi/lam
l = 1/k

radius1, radius2 = args.R1, args.R2

x1, y1, z1 = args.x1, args.y1, args.z1
position1 = [x1/l, y1/l, z1/l]

x2, y2, z2 = args.x2, args.y2, args.z2
position2 = [x2/l, y2/l, z2/l]


beam = Beam()
beam.set_CB(CB)
beam.set_Stokes(stokes)  # set Stokes parameter for beam

bead_list = []

bead = Particle('bead1')
bead.set_radius(radius1/l)  # radius / 1500 nm 
bead.set_n(n1)
bead.set_position(position1)
bead_list.append(bead)

bead2 = Particle('bead2')
bead2.set_radius(radius2/l)   # radius / 10 nm
bead2.set_n(n2)
bead2.set_position(position2)
bead_list.append(bead2)

cluster = MSTM()

cluster.set_beam(beam)  # add beam 
cluster.set_parameter('medium_real_ref_index', n.real)
cluster.set_parameter('medium_imag_ref_index', n.imag)

for bead in bead_list:
    cluster.add_particle(bead)
