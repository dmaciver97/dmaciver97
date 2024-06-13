import numpy as np
import matplotlib.pyplot as plt

gamma = 79.97e-3
v_mol = 1e-30
kT = 1.38e-23*298
r_list = np.linspace(0, 1.6e-10, 1000)
r_star1 = 2*v_mol*gamma/(kT*np.log(1.5))
r_star2 = 2*v_mol*gamma/(kT*np.log(1.8)) 
def free_energy(r, S):
	surface = 4*np.pi*gamma*r**2
	volume = (4*np.pi*r**3/(3*v_mol))*kT*np.log(S) 
	return (surface - volume)/kT, surface/kT, volume/kT

plt.plot(r_list, [free_energy(r,1.5)[0] for r in r_list],
		 c='blue', label = r'$\Delta G_{tot}$ S = 1.5' )
plt.plot(r_list, [free_energy(r,1.8)[0] for r in r_list],
		 linestyle= '--', c='blue',label =r'$\Delta G_{tot}$, S = 1.8' )
plt.plot(r_list, [free_energy(r,1.5)[1] for r in r_list], 
		 c='orange', label =r'$\Delta G_{inf}$')
plt.plot(r_list, [-free_energy(r,1.5)[2] for r in r_list],
		 linestyle= '--',c='green', label =r'$\Delta G_{vol}$, S = 1.5')
plt.plot(r_list, [-free_energy(r,1.8)[2] for r in r_list], 
		 c= 'green', label =r'$\Delta G_{vol}$, S = 1.8')
plt.vlines(r_star1, 0, free_energy(r_star1,1.5)[0], colors='black')
plt.vlines(r_star2, 0, free_energy(r_star2,1.8)[0], colors='black')
plt.axhline(0, c='black', )
plt.legend()
plt.xlim(0, 1.6e-10)
plt.ylabel(r'$\frac{\Delta G}{k_B T}$', fontsize = 'large')
plt.xlabel('Nucleus radius', fontsize = 'large')
plt.xticks(ticks = [r_star1, r_star2], labels = [r'$r_1^*$', r'$r_2^*$'])
plt.yticks(ticks = [0], labels = [0])
plt.show()
plt.savefig('Free_energy_diagram.png')
plt.close()