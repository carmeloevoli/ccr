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

fig = plt.figure(figsize=(8.1, 7.9))
ax = fig.add_subplot(1, 1, 1)

for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)

ax.minorticks_on()
ax.tick_params('both', length=10, width=1.5, which='major', pad=6)
ax.tick_params('both', length=10, width=0, which='minor', pad=6)

plt.xticks(size=30)
plt.yticks(size=30)
# end plot style options

def plot_ionization():
    plt.ylabel(r'Ionization Rate [Myr$^{-1}$]', size=25)
    plt.ylim([1e-9,1e-2])

#    filename = 'output/test_new_2.0_igm.txt'
#    z, pion, crion = np.loadtxt(filename, skiprows=1, usecols=(0,6,7), unpack=True)
#    plt.plot(z, pion, 'k--')
#    plt.text(7,2e-4,'UV',fontsize=24)
#    plt.plot(z, crion * 1.4, 'r', label=r'$\alpha = 2$')

#    filename = 'output/test_new_2.2_igm.txt'
#    z, pion, crion = np.loadtxt(filename, skiprows=1, usecols=(0,6,7), unpack=True)
#    plt.plot(z, crion * 1.4, 'b', label=r'$\alpha = 2.2$')

#    filename = 'output/test_new_2.5_igm.txt'
#    z, pion, crion = np.loadtxt(filename, skiprows=1, usecols=(0,6,7), unpack=True)
#    plt.plot(z, crion * 1.4, 'g', label=r'$\alpha = 2.5$')

#    plt.legend(loc='upper right',fontsize=24,frameon=False)


def plot_tigm():
    plt.ylabel(r'$T_{\rm IGM}$ [K]', size=25)
#plt.ylim([1e-9,1e-2])

    filename = 'output/test_fin_2.2_igm.txt'
    z, y = np.loadtxt(filename, skiprows=1, usecols=(0,2), unpack=True)
    plt.plot(z, y, 'r')

    filename = 'output/test_new_2.2_igm.txt'
    z, y = np.loadtxt(filename, skiprows=1, usecols=(0,2), unpack=True)
    plt.plot(z, y, 'b--')
    
    plt.plot(z, (1.+z)*2.725, 'k--')

def plot_heating():
    plt.ylabel(r'$\Delta T_{\rm IGM}$ [K]', size=25)
    plt.ylim([1e0,1e4])

    kB = 1.3806488e-16 # kelvin / erg
    H0 = 0.7 * 3.2407e-18 # s-1
    OMm = 0.3175
    OMr = 8.6e-5
    OMl = 1. - OMm
    
    filename = 'output/test_new_2.0_igm.txt'
    z, y = np.loadtxt(filename, skiprows=1, usecols=(0,10), unpack=True)
    Hz = H0 * np.sqrt(OMm * (1. + z)**3 + OMr * (1. + z)**4 + OMl);
    dT = 2. / 3. / kB / Hz * np.array(y) / 3.14e13 # kelvin / erg * s * erg / Myr
    plt.plot(z, 1.2 * dT, 'r', label=r'$\alpha = 2$')
    
    filename = 'output/test_new_2.2_igm.txt'
    z, y = np.loadtxt(filename, skiprows=1, usecols=(0,10), unpack=True)
    Hz = H0 * np.sqrt(OMm * (1. + z)**3 + OMr * (1. + z)**4 + OMl);
    dT = 2. / 3. / kB / Hz * np.array(y) / 3.14e13 # kelvin / erg * s * erg / Myr
    plt.plot(z, 1.2 * dT, 'b', label=r'$\alpha = 2.2$')

#filename = 'output/test_new_2.5_igm.txt'
#    z, y = np.loadtxt(filename, skiprows=1, usecols=(0,10), unpack=True)
#    dT = 2. / 3. / kB / Hz * np.array(y) / 3.14e13 # kelvin / erg * s * erg / Myr
#    plt.plot(z, 1.2 * dT, 'g', label=r'$\alpha = 2.5$')

#   plt.plot(z,(1.+z)*2.725,'k--')
    
    #    plt.plot(z,((1.+z)/49.)**2.*40.,'k:')
    
    #plt.text(7,8,r'$T_{\rm CMB}$',fontsize=23)

#plt.legend(loc='upper right',fontsize=21,frameon=False)


#Ionization
ax1 = plt.subplot(211)
#ax1.minorticks_on()
ax1.tick_params('both', length=10, width=1.6, which='major', pad=6)
ax1.tick_params('both', length=0, width=1.3, which='minor', pad=6)

plt.yscale('log')
#plt.axis()#[1e-3,10,1e-2,1e5],interpolation='none')
plt.xlim([6,15])
plt.setp(ax1.get_yticklabels(), fontsize=26)
plt.setp(ax1.get_xticklabels(), visible=False)

#plot_ionization()

plot_tigm()

#Heating
ax2 = plt.subplot(212, sharex=ax1)
#ax2.minorticks_on()
ax2.tick_params('both', length=10, width=1.6, which='major', pad=6)
ax2.tick_params('both', length=0, width=1.3, which='minor', pad=6)

plt.xlabel(r'$z$', size=30)
plt.yscale('log')
plt.xlim([6,15])
plt.setp(ax2.get_yticklabels(), fontsize=26)
plt.setp(ax2.get_xticklabels(), fontsize=26)

plot_heating()

plt.show()

#plt.savefig('cr_igm_rates.pdf', format='pdf', dpi=300)

#plt.savefig('ionization_rate.pdf', format='pdf', dpi=300)
