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
rcParams['lines.linewidth'] = 3
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

def read_file(datafile,xcol,ycol):
    x = []
    y = []
    f = open(datafile,'r')
    header = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        x.append(float(columns[xcol]))
        y.append(float(columns[ycol]))
    f.close()
    data=[]
    data.append(np.array(x))
    data.append(np.array(y))
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

plt.yscale('log')
plt.xscale('log')

plt.xlabel(r'$E$ [GeV]', size=28)
plt.ylabel(r'$t$ [Gyr]', size=28)

plt.axis()#[1e7,1e10,1e-3,1e2])#[1e-3,10,1e-2,1e5],interpolation='none')

E_k = np.logspace(-3,3,100)

plot_timescale(E_k,30.,0.5,'b')

plot_thubble(30,'b')

plot_timescale(E_k,20.,0.05,'r')

plot_thubble(20,'r')



#plt.legend(loc='upper right')

plt.show()

#plt.savefig('timescales.pdf', format='pdf', bbox_inches='tight', dpi=300)