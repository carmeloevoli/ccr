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
    f = open('output/'+datafile,'r')
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

plt.yscale('log')
plt.xscale('log')

plt.xlabel(r'$E_k$ [GeV]', size=28)
plt.ylabel(r'$t_i / t_H$', size=28)

plt.title(r'$z = 20$ - $n_H = 1.7e-3$ - $n_e = 2e-7$')

plt.axis([1e-3,10,1e-2,1e5],interpolation='none')

t_H = read_file('timescales_at_z20.txt',0,1)
t_I = read_file('timescales_at_z20.txt',0,2)
t_C = read_file('timescales_at_z20.txt',0,3)
t_a = read_file('timescales_at_z20.txt',0,4)
t_f = read_file('timescales_at_z20.txt',0,5)

plt.plot(t_H[0],t_I[1]/t_H[1],'r',label='Ionization')
plt.plot(t_H[0],t_C[1]/t_H[1],'b',label='Coulomb')
plt.plot(t_H[0],t_a[1]/t_H[1],'g',label='Adiabatic')
plt.plot(t_H[0],t_f[1]/t_H[1],'m',label='Fragmentation')

t_H = read_file('timescales_at_z10.txt',0,1)
t_I = read_file('timescales_at_z10.txt',0,2)
t_C = read_file('timescales_at_z10.txt',0,3)
t_a = read_file('timescales_at_z10.txt',0,4)
t_f = read_file('timescales_at_z10.txt',0,5)

plt.plot(t_H[0],t_I[1]/t_H[1],'r--')
plt.plot(t_H[0],t_C[1]/t_H[1],'b--')
plt.plot(t_H[0],t_a[1]/t_H[1],'g--')
plt.plot(t_H[0],t_f[1]/t_H[1],'m--')

plt.legend(loc='upper right')

#plt.show()

plt.savefig('timescales.pdf', format='pdf', bbox_inches='tight', dpi=300)