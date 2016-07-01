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

plt.yscale('log')
#plt.xscale('log')

plt.xlabel(r'$z$', size=28)
plt.ylabel(r'$t_{\rm H} / t_{\rm i}$', size=28)

plt.axis()#[1e-3,10,1e-2,1e5],interpolation='none')

t_H = read_file('output/test_no_CR_losses.txt',0,10)
#plt.plot(t_H[0],t_H[1]/t_H[1],'b:',lw=2)

t_I = read_file('output/test_no_CR_losses.txt',0,2)
plt.plot(t_I[0],t_I[1]/t_H[1],'r:')

t_C = read_file('output/test_no_CR_losses.txt',0,3)
plt.plot(t_C[0],t_C[1]/t_H[1],'r--')

t_a = read_file('output/test_no_CR_losses.txt',0,5)
plt.plot(t_a[0],t_a[1]/t_H[1],'b--')

t_pp = read_file('output/test_no_CR_losses.txt',0,4)

t_tot = 1. / (1. / t_C[1] + 1. / t_I[1])
plt.plot(t_C[0],t_tot/t_H[1],'r',label='$E = 1$ MeV')

t_I = read_file('output/test_no_CR_losses.txt',0,6)
plt.plot(t_I[0],t_I[1]/t_H[1],'g:')

t_C = read_file('output/test_no_CR_losses.txt',0,7)
plt.plot(t_C[0],t_C[1]/t_H[1],'g--')

t_a = read_file('output/test_no_CR_losses.txt',0,9)
#plt.plot(t_a[0],t_a[1]/t_H[1],'g')

t_pp = read_file('output/test_no_CR_losses.txt',0,8)

t_tot = 1. / (1. / t_C[1] + 1. / t_I[1])
plt.plot(t_C[0],t_tot/t_H[1],'g',label='$E = 10$ MeV')

plt.text(16.5,10,'Coulomb',size=20,rotation=65)
plt.text(6.5,10,'Ionization',size=20,rotation=-75)
plt.text(20.1,0.06,'Total',size=20,rotation=-10)
plt.text(23,1.8,'Adiabatic',size=20,rotation=0)

plt.ylim([1e-2,1e2])

plt.legend(loc='upper right')

#plt.show()

plt.savefig('timescales_z.pdf', format='pdf', bbox_inches='tight', dpi=300)