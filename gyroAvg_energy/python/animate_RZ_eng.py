"""
10/15/2020
Wes Johnson
This python script animates the orbits of a particle in magnetic feild
"""
#imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import animation
import time 
          
start= time.perf_counter()
#input data
data = np.loadtxt("../data/guidingCenter.txt");times,vx,vy,vz,\
        px,py,pz,eng=data[:,0],data[:,1],data[:,2],data[:,3],\
        data[:,4],data[:,5],data[:,6],data[:,7]

data = np.loadtxt("../../gyroAvg/data/lorenzIon_subRad.txt");timesB,vxB,vyB,vzB,\
        pxB,pyB,pzB,engBor=data[:,0],data[:,1],data[:,2],data[:,3],\
        data[:,4],data[:,5],data[:,6],data[:,7]

end= time.perf_counter()
print(f"Time to Load: {end-start:0.8f}[S]")



#define coordiante transformation: 
def coordsCylTor(px,py,pz):
    R0 = 8.0
    rho = np.sqrt(px*px+py*py)
    phi = np.arctan2(py,px)
    theta = np.arctan2(pz,rho-R0)
    r = np.sqrt((rho-R0)*(rho-R0)+pz*pz)
    return rho,phi,theta,r

"""
#transform coordinates: 
rho,phi,theta,r = coordsCylTor(px,py,pz)
rhoB,phiB,thetaB,rB = coordsCylTor(pxB,pyB,pzB)

# use LaTeX fonts in the plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure()
gs = fig.add_gridspec(1, 2)
ax = fig.add_subplot(gs[0, 0])
axx = fig.add_subplot(gs[0, 1],projection='3d')

ax.set_title(r"\textbf{z vs. $\rho$}")
ax.plot(rho,pz,label = 'Guiding Center',zorder=2)
ax.plot(rhoB,pzB,label = 'Lorenz Ion',zorder=1)
ax.set_ylabel('z [m]')
ax.set_aspect('equal')
ax.set_xlabel(r'$\rho$ [m]')
ax.legend()

axx.plot(px,py,pz)
#make the torus:
n = 100

thetaTor = np.linspace(0, 2.*np.pi, n)
phiTor = np.linspace(0, 2.*np.pi, n)
thetaTor, phiTor = np.meshgrid(thetaTor, phiTor)
c, a = 2000, 500
 
#draw a torus
x = (c + a*np.cos(thetaTor)) * np.cos(phiTor)
y = (c + a*np.cos(thetaTor)) * np.sin(phiTor)
z = a * np.sin(thetaTor)

axx.set_zlim(-c,c)
axx.plot_surface(x, y, z, rstride=5, cstride=5, color='k', edgecolors='w',\
	alpha=.1)
axx.view_init(36, 26)

plt.show()
"""

#do the animation:
FRAMES = 200

#transform coordinates: 
rho,phi,theta,r = coordsCylTor(px,py,pz)
rhoB,phiB,thetaB,rB = coordsCylTor(pxB,pyB,pzB)

positions   = np.array([rho,pz])
positions1  = np.array([rhoB,pzB])
timeSim     = times
q           = 1.60217*10**-19#electron charge 
energySim   = eng/q#units of eV
energyBSim  = engBor/q


#here we animate the resulting trajectories:

fig = plt.figure()
ax = fig.add_subplot()

# use LaTeX fonts in the plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
 
def init():
    ax.set_title(r"\textbf{z vs. $\rho$}")
    ax.set_ylabel('z [m]')
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\rho$ [m]')

def animate(i):
    current_index = int(positions.shape[1] / FRAMES * i)
    ax.cla()
    
    #animate the first trajectory:
    ax.plot(positions[0,:current_index],
              positions[1,:current_index],color='r',zorder=2,
              alpha = 0.7)
    ax.plot([positions[0,current_index]],
            [positions[1,current_index]],
            markerfacecolor='r', markeredgecolor='k',
            marker='o', markersize=5, alpha=1.0,zorder=4)
    #animate the 2nd trajectory:
    ax.plot(positions1[0,:current_index],
              positions1[1,:current_index],color='g',zorder=1,
              alpha = 0.5)
    ax.plot([positions1[0,current_index]],
            [positions1[1,current_index]],
            markerfacecolor='g', markeredgecolor='k',
            marker='o', markersize=5, alpha=1.0,zorder=3)
    ax.set_title(r"\textbf{z vs. $\rho$} with unequal timesteps")
    ax.set_ylabel('z [m]')
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\rho$ [m]')
    ax.annotate(f'Time: {timeSim[current_index]*10**6:0.2f}'+r'[$\mu$S]'+'\n' 
#                 +r'$\delta E_{G.C.}$: '+f'{energySim[0]-energySim[current_index]:0.4f}[eV]'
#                 +'\n'+r'$\delta E_{B.M.}$: '+f'{energyBSim[0]-energyBSim[current_index]:0.4f}[eV]'
#                 +'\n'+r'$E_{G.C.}-E_{B.M.}$: '+
#                 f'{energySim[current_index]-energyBSim[current_index]:0.4f}[eV]'+'\n'+
#                 r'$E_{B.M.}$: '+f'{energyBSim[current_index]/1000:0.4f}[keV]',
                ,xy=(0.7, 0.5), xycoords='axes fraction', fontsize=8,
                horizontalalignment='right', verticalalignment='bottom')
    if i % 50 == 0:
        print(f'{i/FRAMES*100:0.2f}','%')
# call the animator.
print("")
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=FRAMES, interval=100)
#we can save the animation if we uncomment the code below:
# Set up formatting for the movie files
start= time.perf_counter()
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=3600)
anim.save('../plots/GCandBM_RZplane_gyroAvg_eng.mp4', writer=writer)
end= time.perf_counter()
print(f"Time to Animate and Save: {end-start:0.3f}[S]")
