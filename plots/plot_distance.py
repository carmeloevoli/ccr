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


plt.yscale('log')
plt.xscale('log')

plt.xlabel(r'$M$ [M$_\odot$]', size=28)
plt.ylabel(r'$d$ [Mpc]', size=28)

plt.axis([1e7,1e10,1e-3,1e2])#[1e-3,10,1e-2,1e5],interpolation='none')

data = read_file('hmf_20.txt',1,2)

plt.plot(data[0],3e-2*(data[0]*data[1])**(-1./3.)*(1.+20)/21.,'r',label='z=20')

plot_minmass(20.,'r')

data = read_file('hmf_30.txt',1,2)

plt.plot(data[0],3e-2*(data[0]*data[1])**(-1./3.)*(1.+30)/21.,'b',label='z=30')

plot_minmass(30.,'b')

plt.legend(loc='upper right')

#plt.show()

plt.savefig('distance.pdf', format='pdf', bbox_inches='tight', dpi=300)