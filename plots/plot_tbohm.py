#!/bin/bash/python
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import numpy as np

# begin plot style options
rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Helvetica Neue')
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)
rcParams['legend.numpoints'] = 1
rcParams['lines.linewidth'] = 4
rcParams['figure.autolayout'] = True

fig = plt.figure(figsize=(8.1, 7.8))
ax = fig.add_subplot(1, 1, 1)

for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)

ax.minorticks_on()
ax.tick_params('both', length=15, width=1.5, which='major', pad=6)
ax.tick_params('both', length=10, width=1.3, which='minor', pad=6)

plt.xticks(size=28)
plt.yticks(size=28)
# end plot style options

def read_file(datafile,xcol,ycol,zcol):
    x = []
    y = []
    z = []
    f = open(datafile,'r')
    header = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        x.append(float(columns[xcol]))
        y.append(float(columns[ycol]))
        z.append(float(columns[zcol]))
    f.close()
    data=[]
    data.append(np.array(x))
    data.append(np.array(y))
    data.append(np.array(z))
    return data

def plot_minmass(z,color):
    m = 1e8 * (10. / (1. + z))**(1.5)
    plt.plot([m,m],[1e-10,1e10],linestyle=':',color=color)

def plot_timescale(E_k,z,d,color):
    m_p = 0.984
    E_t = E_k + m_p
    p = np.sqrt(E_t**2. - m_p**2.)
    
    DB = 1.1 * p * (1. + z)**(-2.)
    t = d * d / DB
    plt.plot(E_k,t,color=color)

def plot_thubble(z,color):
    t = 19. * (1. + z)**(-3./2.)
    plt.plot([1e-5,1e5],[t,t],linestyle=':',color=color)

def oplot_d_halo():
    data = read_file("hmf_min.txt",0,1,2)
    plt.plot(data[0],3e-2*(data[1]*data[2])**(-1./3.)*(1.+20)/21.,color='b',linestyle='--',label=r'$\langle d \rangle$ halos')

    plt.text(4,0.008,'M$_{min}$',fontsize=25,color='b')

def oplot_CR():
    z = np.linspace(0,30,100)


def oplot_X_ray():
    z = np.linspace(0,30,100)
    
    E = 100.
    d = 0.1 * ((1. + z) / 21.)**(-3) * (E / 300.)**(3.2)
    plt.plot(z,d,'b:',label='X-ray')

    E = 350.
    d = 0.1 * ((1. + z) / 21.)**(-3) * (E / 300.)**(3.2)
    plt.plot(z,d,'b:')

    plt.text(5,10,'$300$ eV',fontsize=24,color='b')
    plt.text(21,0.003,'$100$ eV',fontsize=24,color='b')

plt.yscale('log')
#plt.xscale('log')

plt.xlabel(r'$z$', size=28)
plt.ylabel(r'$\lambda$ [Mpc]', size=28)

plt.axis()#[1e7,1e10,1e-3,1e2])#[1e-3,10,1e-2,1e5],interpolation='none')

#data = read_file("output/test_with_CR_1_MeV_losses.txt",0,6,7)

#z = data[0]
#d = np.sqrt(data[2]*data[1]) / 1e3 # Mpc

#plt.plot(z,d,color='r',label='E = 1 MeV')

#data = read_file("output/test_with_CR_1_MeV_losses.txt",0,12,13)

#z = data[0]
#d = np.sqrt(data[2]*data[1]) / 1e3 # Mpc

#plt.plot(z,d,color='g',label='E = 10 MeV')

#oplot_d_halo()

oplot_X_ray()

plt.ylim([1e-3,1e2])

plt.legend(loc='upper right',fontsize=24)

plt.show()

#plt.savefig('distance_bohm.pdf', format='pdf', dpi=300)
