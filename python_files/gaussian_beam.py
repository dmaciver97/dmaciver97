import numpy as np
import scipy.special
from scipy.special import riccati_jn, jvp, sph_harm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from subprocess import check_call

epislon_0 = 8.8541e-12
mu_0 = np.pi*4e-7
lam = 1064e-9
omega = 3e8/lam
k = omega *np.sqrt(epislon_0*mu_0)

theta_list = np.linspace(0, 2*np.pi, 20)
phi_list = np.linspace(0, np.pi, 20)

for frame in range(50):
    break
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel(r'Theta $[^{\circ}]$')
    ax.set_ylabel(r'Phi $[^{\circ}]$')
    ax.set_zlabel(r'$Y_{lm}$')
    print(f'fig{frame:03d}')
    for theta in theta_list:
        for phi in phi_list:
            sum = 0 
            for n in range(frame):
                for m in range(n):
                    sum += sph_harm(m, n, theta, phi)
            ax.scatter(theta, phi, sum, s = 2, c ='blue')
    plt.savefig(f'./tmp/junk_{frame:03d}.png')
    

cmd = ['ffmpeg', '-r', '60', '-f', 'image2', '-s', '1920x1080', '-i', '.\tmp\junk_%03d.png', '-vcodec', 'libx264', '-crf', '25', '-pix_fmt', 'yuv420p', 'test.mp4']
check_call(cmd)



        