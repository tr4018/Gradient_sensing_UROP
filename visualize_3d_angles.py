# visualization of random_3d_rotation

import numpy as np
import random_3d_rotation
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec

data_number = 10000
original_phi = np.zeros(data_number)
original_theta = np.zeros(data_number)
new_phi = np.zeros(data_number)
new_theta = np.zeros(data_number)
new_x = np.zeros(data_number)
new_y = np.zeros(data_number)
new_z = np.zeros(data_number)

for i in range(data_number):
    new_theta[i],new_phi[i] = random_3d_rotation.random_3d_rotation(original_theta[i],original_phi[i])
    new_x[i] = np.cos(new_phi[i])*np.sin(new_theta[i])
    new_y[i] = np.sin(new_phi[i])*np.sin(new_theta[i])
    new_z[i] = np.cos(new_theta[i])

theta_line=np.arange(0,np.pi/2,0.1)
phi_line=5.7*np.ones(len(theta_line))
new_theta_line,new_phi_line = random_3d_rotation.random_3d_rotation(theta_line,phi_line)
x_line = np.cos(phi_line)*np.sin(theta_line)
y_line = np.sin(phi_line)*np.sin(theta_line)
z_line = np.cos(theta_line)
new_x_line = np.cos(new_phi_line)*np.sin(new_theta_line)
new_y_line = np.sin(new_phi_line)*np.sin(new_theta_line)
new_z_line = np.cos(new_theta_line)


#surface_phi, surface_theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
#surface_radius = 0.9
#surface_x = surface_radius*np.cos(surface_theta)*np.sin(surface_phi)
#surface_y = surface_radius*np.sin(surface_theta)*np.sin(surface_phi)
#surface_z = surface_radius*np.cos(surface_phi)

bins = 60
sample_every = 20
fig = plt.figure(figsize=[6.4,6])
spec = gridspec.GridSpec(ncols=2, nrows=2, hspace=0.4, height_ratios=[1,5])
ax1 = fig.add_subplot(spec[0])
ax1.set_xlim([0,np.pi])
ax1.set_xlabel(r'$\theta$')
ax1.set_ylabel(r'occurrence')
ax1.set_xticks([0,np.pi/2, np.pi])
ax1.set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$'])
ax1.hist(new_theta,bins)
ax2 = fig.add_subplot(spec[1])
ax2.set_xlim([0,2*np.pi])
ax2.set_xlabel(r'$\phi$')
ax2.set_xticks([0,np.pi, 2*np.pi])
ax2.set_xticklabels([r'$0$', r'$\pi$', r'$2\pi$'])
ax2.hist(new_phi,bins)
ax3 = fig.add_subplot(spec[2:], projection='3d')
ax3.set_xticks([-1,0,1])
ax3.set_xlabel(r'$x$')
ax3.set_yticks([-1,0,1])
ax3.set_ylabel(r'$y$')
ax3.set_zticks([-1,0,1])
ax3.set_zlabel(r'$z$')
#ax3.plot_surface(surface_x,surface_y,surface_z,  rstride=1, cstride=1, color='c', alpha=1.0, linewidth=0)
ax3.scatter(new_x[1::sample_every],new_y[1::sample_every],new_z[1::sample_every],color="k",s=20)
ax3.plot3D(x_line,y_line,z_line,'blue')
ax3.plot3D(new_x_line,new_y_line,new_z_line,'red')
plt.savefig("RandomPhiTheta01.png")