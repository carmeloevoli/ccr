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
plt.ylabel(r'$\Lambda_{\rm ion}$ [Myr$^{-1}$]', size=30)

plt.axis()#[1e-3,10,1e-2,1e5],interpolation='none')

filename = 'output/test_with_CR_2.5_losses.txt'
z, pion, crion = .. 

#data = read_file('output/test_with_CR_100_keV_igm.txt',0,6)
#plt.plot(data[0],data[1],'b',label='Ph-Ion')

#data = read_file('output/test_with_CR_100_keV_igm.txt',0,7)
#plt.plot(data[0],data[1],'b--',label='CR-Ion')

#data = read_file('output/test_with_CR_100_keV_igm.txt',0,9)
#plt.plot(data[0],data[1],'r',label='Ph-Heating')

#data = read_file('output/test_with_CR_100_keV_igm.txt',0,10)
#plt.plot(data[0],data[1],'r--',label='CR-Heating')

plt.show()

#plt.savefig('photo_ionization_rate.pdf', format='pdf', bbox_inches='tight', dpi=300)
