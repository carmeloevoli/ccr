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

plt.xticks(size=30)
plt.yticks(size=30)
# end plot style options

plt.yscale('log')
#plt.xscale('log')

plt.xlabel(r'$z$', size=30)
plt.ylabel(r'$t_{\rm i} / t_{\rm H}$', size=30)

plt.axis()#[1e-3,10,1e-2,1e5],interpolation='none')

filename = 'output/test_no_CR_losses.txt'
filename = 'output/test_with_CR_2.5_losses.txt'

z, t_I, t_C, t_pp, t_a, t_H = np.loadtxt(filename, skiprows=1, usecols=(0,2,3,4,5,14), unpack=True)

plt.plot(z, t_I / t_H, 'r:')
plt.plot(z, t_C / t_H, 'r--')
plt.plot(z, t_a / t_H, 'b--')

t_tot = 1. / (1. / t_C + 1. / t_I)
plt.plot(z, t_tot / t_H, 'r', label='$E = 1$ MeV')

z, t_I, t_C, t_pp, t_a, t_H = np.loadtxt(filename, skiprows=1, usecols=(0,8,9,10,11,14), unpack=True)

plt.plot(z, t_I / t_H, 'g:')
plt.plot(z, t_C / t_H, 'g--')
#plt.plot(z, t_a / t_H, 'b--')

t_tot = 1. / (1. / t_C + 1. / t_I)
plt.plot(z, t_tot / t_H, 'g', label='$E = 10$ MeV')

plt.text(16.5,10,'Coulomb',size=22,rotation=56)
plt.text(6.5,10,'Ionization',size=22,rotation=-70)
plt.text(15,0.04,'Total',size=22)
plt.text(1,0.9,'Adiabatic',size=22,rotation=0)

plt.ylim([1e-2,1e2])
plt.xlim([0,18])

plt.legend(loc='lower left',fontsize=20,frameon=False)

#plt.show()

plt.savefig('timescales_z.pdf', format='pdf', bbox_inches='tight', dpi=300)
