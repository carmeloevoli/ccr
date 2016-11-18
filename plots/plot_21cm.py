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

plt.xticks(size=32)
plt.yticks(size=32)
# end plot style options

def plot_delta_T():
    #plt.ylabel(r'$\Delta T_{\rm IGM}$ [K]', size=32)
    #plt.ylim([1e0,1e4])

    kB = 1.3806488e-16 # kelvin / erg
    h = 0.7
    H0 = h * 3.2407e-18 # s-1
    OMm = 0.3175
    OMr = 8.6e-5
    OMl = 1. - OMm
    delta = 1.0
    
    filename = 'output/test_new_2.0_igm.txt'
    z, x_II, T_k = np.loadtxt(filename, skiprows=1, usecols=(0,1,2), unpack=True)

    T_k = 1.2 * T_k # Helio

    T_CMB = (1.+z)*2.725 # K
    Hz = H0 * np.sqrt(OMm * (1. + z)**3 + OMr * (1. + z)**4 + OMl);
    #dT = 2. / 3. / kB / Hz * np.array(y) / 3.14e13 # kelvin / erg * s * erg / Myr

    T = (1. / h) * (1. - x_II) * (1. + delta) * np.sqrt(((1. + z) / 10.) * (0.3 / OMm))

    ya_eff = 0.
    T_star = 0.068 # K
    A_10 = 2.85e-15 # s-1

    n_b = 2.54e-7 / (1. + z)**3.0
    n_H = n_b
    n_p = x_II * n_H
    n_e = n_p
    
    k = 3.1e-11 * T_k**0.357 * np.exp(-32./T_k) # cm3/s
    log_gamma_e = -9.607 + 0.5 * np.log10(T_k) * np.exp(-(np.log10(T_k)**4.5/1800.))
    gamma_e = 10.0**log_gamma_e # cm3/s

    C_H = n_H * k
    C_e = n_e * gamma_e
    C_p = n_p * k * 3.2

    y_c = T_star / A_10 / T_k * (C_H + C_e + C_p)
    
    T_s = (T_CMB + (ya_eff + y_c) * T_k) / (1. + ya_eff + y_c)
    
    dT = 0.016 * T * (1. - T_CMB / T_s)

    print min(dT), max(dT)
    
#    plt.plot(z, dT, 'r', label=r'$\alpha = 2$')

    plt.plot(z, T_CMB, label='cmb')
    plt.plot(z, T_k, label='cmb')

    #plt.plot(z,(1.+z)*2.725,'k--')
    #plt.text(13,65,r'$T_{\rm CMB}$',fontsize=32)
    #plt.legend(loc='upper right',fontsize=24,frameon=False)

plt.yscale('log')
plt.xlabel(r'$z$', size=32)
plt.axis()#[1e-3,10,1e-2,1e5],interpolation='none')
plt.xlim([6,15])

plot_delta_T()

plt.show()

#plt.savefig('heating_rate.pdf', format='pdf', dpi=300)


