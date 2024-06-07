#* define box of dimensions given below
x_range = 2e-6
y_range = 2e-6
z_range = 2e-6
from scipy.fft import ifftn, fftn

#* each axis is divided into segments of length i_range/Ni
#* gridding can be increased for increased accuracy
Nx = 12
Ny = 12
Nz = 12

#* X spacing
dXx = 2*x_range/Nx
dXy = 2*y_range/Ny
dXz = 2*z_range/Nz

#* q_kx = 2pi * kx * Nx^2/length_x
dqx = 1/x_range #! units of m^-1
dqy = 1/y_range #! units of m^-1
dqz = 1/z_range #! units of m^-1
#? dj = scan_range/N

x_axis = np.linspace(-x_range, x_range, Nx)
y_axis = np.linspace(-y_range, y_range, Ny)
z_axis = np.linspace(-z_range, z_range, Nz)

pos_array = np.zeros((Nx, Ny, Nz), dtype=np.ndarray)
force_array = np.zeros((Nx, Ny, Nz), dtype=np.ndarray)
q_array = np.zeros((Nx, Ny, Nz), dtype=np.ndarray)

#* scan over entire X space and compute x, y, z force components
for jx, x in enumerate(x_axis):
    for jy, y in enumerate(y_axis):
        for jz, z in enumerate(z_axis):
            #* assign position and force components according to j index
            position = [x, y, z]
            pos_array[jx, jy, jz] = position
           
            force = force_func(pos=position,
                               rot_matrix=as_rotation_matrix(q))
            force_array[jx, jy, jz] = [force[0], force[1], force[2]]

#fig = plt.figure()
#ax = fig.add_subplot(projection='3d')
#for i, Ax1 in enumerate(force_array):
#    for j, Ax2 in enumerate(Ax1):
#        for k, Fq in enumerate(Ax2):      
#            pos = pos_array[i,j,k]
#            x = pos[0]*1e6
#            y = pos[1]*1e6
#            z = pos[2]*1e6
#            f_grad = force_array[i,j,k]
#            if f_grad != 0:
#                fx = f_grad[0]
#                fy = f_grad[1]
#                fz = f_grad[2]
#                ax.quiver(x, y, z, fx,fy,fz, length = 0.5, normalize = True)
#            else:
#                pass
#
#ax.set_xlim(-x_range*1e6, x_range*1e6)
#ax.set_ylim(-y_range*1e6, y_range*1e6)
#ax.set_zlim(-z_range*1e6, z_range*1e6)
#
#plt.show()
#plt.close()

force_qx = np.zeros((Nx,Ny,Nz))
force_qy = np.zeros((Nx,Ny,Nz))
force_qz = np.zeros((Nx,Ny,Nz))

force_array_q = np.zeros((Nx,Ny,Nz), dtype=np.ndarray)

#* Assign point in q space using k indexing
for kx  in range(Nx):
    for ky in range(Ny):
        for kz in range(Nz):
            #* determine q vector at that particular point
            #* assign by q vector and inverse transformed force using same k index
            q_vec = np.array([dqx*kx, dqy*ky, dqz*kz])
            q_array[kx, ky, kz] = q_vec
            Fqx = 0
            Fqy = 0
            Fqz = 0
            for jx, Ax1 in enumerate(force_array):
                for jy, Ax2 in enumerate(Ax1):
                    for jz, Force in enumerate(Ax2):
                        #* sum over entire force field at point q_kx, q_ky, q_kz

                        fx = Force[0]
                        fy = Force[1]
                        fz = Force[2]
                        pos = pos_array[jx, jy, jz]

                        dot_x = q_vec[0]*pos[0]
                        dot_y = q_vec[1]*pos[1]
                        dot_z = q_vec[2]*pos[2]

                        Fqx += fx * np.exp(-1j*dot_x-1j*dot_y-1j*dot_z)
                        Fqy += fy * np.exp(-1j*dot_x-1j*dot_y-1j*dot_z)
                        Fqz += fz * np.exp(-1j*dot_x-1j*dot_y-1j*dot_z)

            #* Multiply by physical spacing and divide by (2pi)^0.5
            Fqx *= dXx / (2*np.pi)**0.5 * (dXy / (2*np.pi)**0.5) * (dXz / (2*np.pi)**0.5)
            Fqy *= dXx / (2*np.pi)**0.5 * (dXy / (2*np.pi)**0.5) * (dXz / (2*np.pi)**0.5)
            Fqz *= dXx / (2*np.pi)**0.5 * (dXy / (2*np.pi)**0.5) * (dXz / (2*np.pi)**0.5)

            force_array_q[kx, ky, kz] = [Fqx, Fqy, Fqz]

print('Inverting Fourier Transform')
def get_grad(q, Fq):
    q2 = np.dot(q,q)
    if Fq == 0 or q2 == 0:
        return [0,0,0]
    else:
        intergrand = [q[0]*np.dot(q, Fq)/(q2*(2*np.pi)**1.5),
                      q[1]*np.dot(q, Fq)/(q2*(2*np.pi)**1.5),
                      q[2]*np.dot(q, Fq)/(q2*(2*np.pi)**1.5)]
        return intergrand  

def get_scat(q, Fq):
    q2 = np.dot(q,q)
    if Fq == 0 or q2 == 0:
        return [0,0,0]
    else:
        intergrand = np.cross(np.cross(q, Fq), q)/(q2*(2*np.pi)**1.5)
        return intergrand
   
GF_array = np.zeros((Nx,Ny,Nz), dtype=np.ndarray)
SF_array = np.zeros((Nx,Ny,Nz), dtype=np.ndarray)

#* like before we now use the j indexing system to define our position vector
#* this is used to also index our gradient and scattering forces
for jx in range(Nx):
    for jy in range(Ny):
        for jz in range(Nz):
            pos = pos_array[jx, jy, jz]

            F_grad = [0,0,0]
            F_scat = [0,0,0]
           
            for kx, Ax1 in enumerate(force_array_q):
                for ky, Ax2 in enumerate(Ax1):
                    for kz, Fq in enumerate(Ax2):
                        q_vec = q_array[kx,ky,kz]
                       
                        q2 = np.hypot(q_vec[0], np.hypot(q_vec[1], q_vec[2]))**2 #! magnitude of q^2
                        q2 = np.dot(q_vec, q_vec)
                        dot_x = q_vec[0]*pos[0]
                        dot_y = q_vec[1]*pos[1]
                        dot_z = q_vec[2]*pos[2]

                        grad_int = get_grad(q_vec, q2, Fq)
                        scat_int = get_scat(q_vec, q2, Fq)

                        F_grad[0] += grad_int[0]*np.exp(1j*dot_x+1j*dot_y + 1j*dot_z)
                        F_grad[1] += grad_int[1]*np.exp(1j*dot_x+1j*dot_y + 1j*dot_z)
                        F_grad[2] += grad_int[2]*np.exp(1j*dot_x+1j*dot_y + 1j*dot_z)

                        F_scat[0] += scat_int[0]*np.exp(1j*dot_x+1j*dot_y + 1j*dot_z)
                        F_scat[1] += scat_int[1]*np.exp(1j*dot_x+1j*dot_y + 1j*dot_z)
                        F_scat[2] += scat_int[2]*np.exp(1j*dot_x+1j*dot_y + 1j*dot_z)

            F_grad[0] *= dqx*dqy*dqz
            F_grad[1] *= dqx*dqy*dqz
            F_grad[2] *= dqx*dqy*dqz

            F_scat[0] *= dqx*dqy*dqz
            F_scat[1] *= dqx*dqy*dqz
            F_scat[2] *= dqx*dqy*dqz
           
            GF_array[jx, jy, jz] = [F_grad[0].real/(Nx*Ny*Nz), F_grad[1].real/(Nx*Ny*Nz), F_grad[2].real/(Nx*Ny*Nz)]
            SF_array[jx, jy, jz] = [F_scat[0].real/(Nx*Ny*Nz), F_scat[1].real/(Nx*Ny*Nz), F_scat[2].real/(Nx*Ny*Nz)]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for i, Ax1 in enumerate(GF_array):
    for j, Ax2 in enumerate(Ax1):
        for k, Fq in enumerate(Ax2):      
            pos = pos_array[i,j,k]
            x = pos[0]*1e6
            y = pos[1]*1e6
            z = pos[2]*1e6
            f_grad = GF_array[i,j,k]
            if f_grad != 0:
                fx = f_grad[0]
                fy = f_grad[1]
                fz = f_grad[2]
                ax.quiver(x, y, z, fx,fy,fz, length = 0.25, normalize = True)
            else:
                pass

ax.set_xlim(-x_range*1e6, x_range*1e6)
ax.set_ylim(-y_range*1e6, y_range*1e6)
ax.set_zlim(-z_range*1e6, z_range*1e6)

#ax.set_zlim(-z_range*1e6, z_range*1e6)

plt.show()
