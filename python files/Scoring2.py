#!/usr/bin/python3.9

import numpy as np
from statistics import mean, stdev
from scipy.interpolate import splrep, splev, splint
import quaternion
from subprocess import check_call
import random as rand
from quaternion import as_rotation_matrix
from multiprocessing import Process, Queue, Lock
import matplotlib.pyplot as plt
import Experiment3
import xml.etree.ElementTree as ET
from datetime import datetime
 
class Scoring2(object):

    def __init__(self):
        self.traj = np.load('traj_test.csv.npy')
        
        self.time_steps = 1000
        self.dt = 1e-3
        self.pkl = 'fit.pkl'
        self.sim = Experiment3.Experiment3(self.pkl, angle_list=None)
        self.angle_list = [15,55,90]

    def test(self, angle_list, epsilon):
        #Main function for iterating through the main traj list and assigning an estimtated orienation
        self.angle_list = angle_list
        self.sim.set_measurement_angles(self.angle_list)
        self.t_list, self.dist_list, self.est_list = [], [], []
        self.prob_max_list =[]
        s_list, est, true=[], [], []
        #print(self.sim.n_list)
        

        for count, frame in enumerate(self.traj[0:10]):
            #print(count)
            q = np.quaternion(frame[3], frame[4], frame[5], frame[6])
            R = as_rotation_matrix(q)
            s = np.array([R[0,2], R[1,2], R[2,2]])
            s_list.append(s)
            
            fit = self.sim.get_scatter_fit(s)
            I_raw = [splev(angle, fit) for angle in self.angle_list]

            index  = self.sim.get_distances(s)
            n = self.sim.n_list[index]['n']
            true.append(n)
            #s_prior is used to inform our prior estimate of the particle's orienation
            if count != 0:
                s_prior = n_est
                
            else:
                s_prior = s

            index_est, self.prob_list = self.sim.estimate_orientation(I_raw,
                                                                      s_prior,
                                                                      epsilon)
        
            n_est =  self.sim.n_list[index_est]['n']
            est.append(n_est)

            self.prob_max_list.append(self.prob_list[index])
            #if n[0] == n_est[0] and n[1] == n_est[1] and n[2] == n_est[2]:
             #   self.prob_max_list.append(self.prob_list[index])
            #else:
             #  self.prob_max_list.append(self.prob_list[index_est]-
                                         #self.prob_list[index])

            dist = self.sim.get_dist(s, n)
            dist_est = self.sim.get_dist(s, n_est)
            #print(s)
            #print(n, dist)
            #print(n_est, dist_est)
            self.dist_list.append(dist)
            self.est_list.append(dist_est)
            self.t_list.append(self.dt * count)

        plt.hist(self.dist_list, bins =200, density =True)
        plt.show()

        np.save('true.csv', true)
        np.save('est.csv', est)

        now = datetime.now()
        current_time = now.strftime('%H_%M_%S_%d%m%Y')
        np.save(f'dump/dist{self.angle_list}{epsilon}.csv', self.dist_list)
        np.save(f'dump/est{self.angle_list}{epsilon}.csv', self.est_list)
        return self.score_run()#, self.dist_list, self.est_list

    def score_run(self):
        #Calculate and return the divergence of a orientation estimation
        K_l = 0

        for x in range(len(self.prob_max_list)):
            if self.prob_max_list[x] == 0:
                K_l += 5
            else:
                K_l += np.log(1/self.prob_max_list[x])

        worst = len(self.prob_max_list)*(np.log(30))

        self.improve = worst/K_l

        print(self.improve)

        return self.improve

    def plot_prob_dist(self, index, index_est):
        plt.plot(range(1,31), self.prob_list)
        plt.scatter(index_est+1, self.prob_list[index_est],
                    label = 'Our estimate')
        plt.scatter(index+1, self.prob_list[index],
                    label ='Correct Result')
        plt.legend(fontsize =10, loc ='lower right')
        plt.xlabel('Reference Oreientation no.')
        plt.ylabel('Probability ($p(y_{exp}|\hat{n}_{\u03B1})$)')
        plt.show()
    
    def show_divergence(self, dist, est, cap, F):
        #provide a visual plot of the divergence between our estimate and best estimate
        #print(self.t_list, self.dist_list)
        dist_list = np.load(dist)
        est_list = np.load(est)
        t_list = [t*1e-3 for t in range(len(dist_list))]
        plt.close('all')
        lab = ['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{3}$',
               r'$\frac{\pi}{2}$', r'$\frac{2}{3}\  \pi$',
               r'$\frac{3}{4}\  \pi$',r'$\pi$']
        loc = [0,0.25*np.pi,1/3*np.pi,
               0.5*np.pi,2/3*np.pi,0.75*np.pi,np.pi]
        now = datetime.now()
        current_time = now.strftime('%H_%M_%S_%d%m%Y')
        plt.plot(t_list, dist_list,
                 label= 'Best Possible Result', zorder=10)
        plt.plot(t_list, est_list,
                 label = 'Estimation Result', zorder=0)
        plt.axhline(y=0.895607, ls = 'dotted')
        plt.title(label = f'Improvement = {round(F, 4)}',
                  loc = 'left', fontsize ='xx-large')
        plt.ylim(0, np.pi)
        plt.yticks(loc, lab, fontsize='x-large')
        plt.ylabel(r'Angle between $\hat{n}$ & $\hat{s}$ [radians]',
                   fontsize='x-large')
        plt.xlabel('Time [s]', fontsize='x-large')
        plt.legend(fontsize = 'xx-large')
        #plt.savefig(f'PhD_Notes/Divergence_disturbutions/Div_dist{current_time}.png', dpi =1000)
        plt.show()

    def create_y_space_labels(self):
        y_dat = []
        for n_dict in self.sim.n_list:
            y_dat.append(n_dict['y'])
        y_dict = {}
        for label, signal in enumerate(y_dat):
            y_dict[str(signal)]=[]
        for label, signal in enumerate(y_dat):
            label = label+1
            y_dict[str(signal)].append(label)

        return y_dict

    def annotated_plot(self, animate):
        if animate == True:
            for i in range(0,200):
                #self.angle_list = [15,45,90]
                #self.sim.n_list = self.sim.set_measurement_angles(self.angle_list)
                y_dict = self.create_y_space_labels()
                print(1.8*i)
                fig1 = plt.figure()
                ax1 = fig1.add_subplot(projection ='3d')
                for label, n_dict in enumerate(self.sim.n_list):
                    n = n_dict['n']
                    ax1.view_init(elev = 10, azim = 1.8*i)
                    ax1.scatter(n[0], n[1], n[2], alpha=1, color ='blue')
                    ax1.text(n[0], n[1], n[2], str(label+1))
                ax1.set_xlabel('$n_x$')
                ax1.set_ylabel('$n_y$')
                ax1.set_zlabel('$n_z$')
                ax1.set_title(label = 'Vector space', fontsize =15,
                              loc='center')
                plt.savefig(f'./PhD_Notes/Plots/RCE files/tmp vec/Vsp{i:03d}.png')
                plt.close()

                fig2 = plt.figure()
                ax2 = fig2.add_subplot(projection ='3d')
                for label, n_dict in enumerate(self.sim.n_list):
                    y = n_dict['y']
                    signal  = [y[0], y[1], y[2]]
                    ax2.view_init(elev = 10, azim = 1.8*i)
                    ax2.scatter(y[0], y[1], y[2], alpha=1, color ='blue')
                    ax2.text(y[0], y[1], y[2], y_dict[str(signal)])
                ax2.set_xlabel(r'$y_{\theta_1}$')
                ax2.set_ylabel(r'$y_{\theta_2}$')
                ax2.set_zlabel(r'$y_{\theta_3}$')
                ax2.set_title(label = 'Intensity space', fontsize =15,
                              loc='center')
                plt.savefig(f'./PhD_Notes/Plots/RCE files/tmp int/Isp{i:03d}.png')
                plt.close()

            cmd1 = ['ffmpeg', '-r', '25', '-f', 'image2', '-s', '1920x1080', '-i',
                   './PhD_Notes/Plots/RCE files/tmp vec/Vsp%3d.png', '-vcodec',
                    'libx264', '-crf', '25', '-pix_fmt', 'yuv420p', 'Vsp.mp4']
            check_call(cmd1)

            cmd2 = ['ffmpeg', '-r', '25', '-f', 'image2', '-s', '1920x1080', '-i',
                   './PhD_Notes/Plots/RCE files/tmp int/Isp%3d.png', '-vcodec',
                    'libx264', '-crf', '25', '-pix_fmt', 'yuv420p', 'Isp.mp4']
            check_call(cmd2)
            

        else:
            kwargs = {'color':'r'}
            y_dict = self.create_y_space_labels()
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(projection ='3d')
            
            for label, n_dict in enumerate(self.sim.n_list):
                n = n_dict['n']
                ax1.scatter(n[0], n[1], n[2], s = 20, alpha=1, color ='blue')
                ax1.text(n[0], n[1], n[2], str(label+1))
            ax1.set_xlabel('$n_x$')
            ax1.set_ylabel('$n_y$')
            ax1.set_zlabel('$n_z$')
            ax1.set_title(label = 'Vector space', fontsize =15, loc='center')
            plt.show()

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(projection ='3d')
            for label, n_dict in enumerate(self.sim.n_list):
                y = n_dict['y']
                signal  = [y[0], y[1], y[2]]
                ax2.scatter(y[0], y[1], y[2], alpha=1, color ='blue')
                ax2.text(y[0], y[1], y[2], y_dict[str(signal)])
            ax2.set_xlabel('$y_x$')
            ax2.set_ylabel('$y_y$')
            ax2.set_zlabel('$y_z$')
            ax2.set_title(label = 'Intensity space', fontsize =15, loc='center')
            plt.savefig('./PhD_Notes/Plots/RCE files/Intensity space')
            
    def optimize(self,i):
        #optimize uses shorter lists of orientations to quickly find a close fit for the angle list
        #using stochastic optimization + simulated annealing techniques. 
        s_list = [self.traj[x] for x in range(i*(100), (i+1)*100)]
        self.traj = s_list
        self.time_steps = len(s_list)

        dtheta = [5,10,11.6]
        theta_max  = [25, 60, 100]
        theta_min  = [10, 30, 65]
        K_l_old = self.K_l
        beta = 1
        fail = 0
        while beta <= 100 :
            print('----')
            angles_old = self.angle_list
            
            angles_new = [angle+rand.uniform(-1,1)*dtheta[x] for angle, x in zip(self.angle_list, range(0,3))]
            
            for x in range(0,3):
                if angles_new[x] >= theta_max[x] or angles_new[x] <= theta_min[x]:
                    angles_new[x] = angles_old[x]

            self.angle_list = angles_new
            
            K_l_new = self.test(angles_new)
            print(K_l_old)

            if rand.random() < np.exp(-beta*(K_l_old-K_l_new)):
                print("Sucess")
                K_l_old = K_l_new
                dtheta = [x*0.5 for x in dtheta]
                beta = beta*1.5
                print(beta, angles_new)
                self.K_l = K_l_new
                fail = 0 
            else:
                fail += 1
                print("Failure")
                self.angle_list = angles_old

            if fail >= 20:
                beta = beta * 1.5
                print(beta)
                fail = 0 
                dtheta = [x*0.25 for x in dtheta]

        print(self.K_l, angles_new)

        return angles_new

    def plot_scatter_ref(self, xml, n):
        xml_file = ET.parse(xml)
        sim = xml_file.getroot()
        for data in sim:
            s_vec = [float(data.attrib['s0']), float(data.attrib['s1']),float(data.attrib['s2'])]
            if s_vec[0]  == n[0] and s_vec[1] == n[1] and s_vec[2] == n[2]:
                #print(s_vec)
                theta_list =[]
                I_raw =[]
                for lines in data:
                    theta_list.append(float(lines.attrib['theta']))
                    I_raw.append(float(lines.attrib['I']))

        lab = ['0', r'$\theta_1$', '20', '40', r'$\theta_2$', '60', '80'
               , r'$\theta_3$', '100', '120', '140', '160', '180']
        loc = [0, self.angle_list[0], 20, 40, self.angle_list[1], 60, 80,
               self.angle_list[2], 100, 120, 140, 160, 180]
        
        plt.plot(theta_list, I_raw)
        for angle in self.angle_list:
            plt.vlines(angle, 1e-5, 1e2, color='black', ls='dotted')
        plt.xlabel(r'$\theta$', fontsize = 15)
        plt.ylabel(r'$I(\theta)$', fontsize = 15)
        plt.yscale('log')
        plt.xlim([0,180])
        plt.ylim([1e-5, 1e2])
        plt.xticks(loc,  lab)
        plt.show()
                


run= Scoring2()
run.test([15,55,90], 0.15)
