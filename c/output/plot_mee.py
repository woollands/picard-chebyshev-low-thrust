import pandas as pd
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.basemap import Basemap
import const as cn
from functions import mee2rv
from functions import drawEarth
from functions import thrust_angle

data=pd.read_csv('data.dat', sep='\t')
t = data.values[:,0];
p = data.values[:,1];
f = data.values[:,2];
g = data.values[:,3];
h = data.values[:,4];
k = data.values[:,5];
L = data.values[:,6];
m = data.values[:,7];
plam = data.values[:,8];
flam = data.values[:,9];
glam = data.values[:,10];
hlam = data.values[:,11];
klam = data.values[:,12];
Llam = data.values[:,13];
mlam = data.values[:,14];

[r,v] = mee2rv(p,f,g,h,k,L,cn.mu)
rho = 1.0;
[u_inert,u_lvlh,S,F,Pa,delta,zeta] = thrust_angle(data.values,rho,True)

fig, (ax1, ax2, ax3) = plt.subplots(3)
fig.suptitle('States 1')
ax1.plot(t,p)
ax1.set_ylabel('p (DU)')
ax2.plot(t,f)
ax2.set_ylabel('f')
ax3.plot(t,g)
ax3.set_xlabel('Time (TU)')
ax3.set_ylabel('g')

fig, (ax1, ax2, ax3) = plt.subplots(3)
fig.suptitle('States 2')
ax1.plot(t,h)
ax1.set_ylabel('h')
ax2.plot(t,k)
ax2.set_ylabel('k')
ax3.plot(t,L)
ax3.set_xlabel('Time (TU)')
ax3.set_ylabel('L')

fig, (ax1, ax2, ax3) = plt.subplots(3)
fig.suptitle('Costates 1')
ax1.plot(t,plam)
ax1.set_ylabel(r'$\lambda_{\mathrm{p}}$')
ax2.plot(t,flam)
ax2.set_ylabel(r'$\lambda_{\mathrm{f}}$')
ax3.plot(t,glam)
ax3.set_xlabel('Time (TU)')
ax3.set_ylabel(r'$\lambda_{\mathrm{g}}$')

fig, (ax1, ax2, ax3) = plt.subplots(3)
fig.suptitle('Costates 2')
ax1.plot(t,hlam)
ax1.set_ylabel(r'$\lambda_{\mathrm{h}}$')
ax2.plot(t,klam)
ax2.set_ylabel(r'$\lambda_{\mathrm{k}}$')
ax3.plot(t,Llam)
ax3.set_xlabel('Time (TU)')
ax3.set_ylabel(r'$\lambda_{\mathrm{L}}$')

fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle('Mass')
ax1.plot(t,m)
ax1.set_ylabel('Mass (kg)')
ax2.plot(t,mlam)
ax2.set_ylabel(r'$\lambda_{\mathrm{m}}$')
ax2.set_xlabel('Time (TU)')
ax1.grid()
ax2.grid()

# fig, (ax1, ax2) = plt.subplots(2)
# fig.suptitle('Power')
# ax1.plot(t,Pa)
# ax1.set_ylabel('Available Power')
# ax2.plot(t,delta)
# ax2.plot(t,zeta)
# ax2.plot(t,delta*zeta)
# ax2.legend([r'$\delta$',r'$\zeta$',r'$\delta$ * $\zeta$'])
# ax2.set_ylabel('Continuation Parameters')
# ax2.set_xlabel('Time (TU)')
# ax1.grid()
# ax2.grid()

fig, (ax1, ax2, ax3) = plt.subplots(3)
fig.suptitle('Power and Thrust')
ax1.plot(t,S)
ax1.set_ylabel('Switch Function')
ax2.plot(t,Pa)
ax2.set_ylabel('Power (W)')
ax3.plot(t,F)
ax3.set_ylabel('Thrust (N)')
ax3.set_xlabel('Time (TU)')
ax1.grid()
ax2.grid()
ax3.grid()

[x, y, z, bm] = drawEarth(1)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(r[:,0],r[:,1],r[:,2])
ax.quiver(r[:,0],r[:,1],r[:,2],u_lvlh[:,0],u_lvlh[:,1],u_lvlh[:,2],length=0.2,colors='red')
ax.plot_surface(x, y, z,  rstride=4, cstride=4, alpha=0.5, facecolors=bm)
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_zlabel('z (km)')
ax.set_xlim(-6,6)
ax.set_ylim(-6,6)
ax.set_zlim(-6,6)
# axis.set_clip_on(self,b)

plt.show()
