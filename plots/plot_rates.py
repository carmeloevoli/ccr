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

def plot_ionization():
    plt.ylabel(r'Ionization Rate [Myr$^{-1}$]', size=28)
    plt.ylim([1e-9,1e-3])

    data = read_file('output/test_with_CR_2.0_igm.txt',0,6)
    plt.plot(data[0],data[1] / 10.,'r--',label='UV')

    data = read_file('output/test_with_CR_2.0_igm.txt',0,7)
    plt.plot(data[0],1.4*data[1],'r',label=r'$\alpha = 2$')

    data = read_file('output/test_with_CR_2.2_igm.txt',0,7)
    plt.plot(data[0],1.4*data[1],'b',label=r'$\alpha = 2.2$')

    data = read_file('output/test_with_CR_2.5_igm.txt',0,7)
    plt.plot(data[0],1.4*data[1],'g',label=r'$\alpha = 2.5$')

    plt.legend(loc='upper right',fontsize=16)

def plot_heating():
    plt.ylabel(r'$\Delta T_{\rm IGM}$ [K]', size=28)
    plt.ylim([1e0,1e4])

    kB = 1.3806488e-16 # kelvin / erg
    H0 = 0.7 * 3.2407e-18 # s-1
    OMm = 0.3175
    OMr = 8.6e-5
    OMl = 1. - OMm
    
    data = read_file('output/test_with_CR_2.0_igm.txt',0,10)
    z = np.array(data[0])
    Hz = H0 * np.sqrt(OMm * (1. + z)**3 + OMr * (1. + z)**4 + OMl);
    dT = 2. / 3. / kB / Hz * np.array(data[1]) / 3.14e13 # kelvin / erg * s * erg / Myr
    plt.plot(data[0],1.4*dT,'r',label=r'$\alpha = 2$')
    
    data = read_file('output/test_with_CR_2.2_igm.txt',0,10)
    dT = 2. / 3. / kB / Hz * np.array(data[1]) / 3.14e13 # kelvin / erg * s * erg / Myr
    plt.plot(data[0],1.4*dT,'b',label=r'$\alpha = 2.2$')
    
    data = read_file('output/test_with_CR_2.5_igm.txt',0,10)
    dT = 2. / 3. / kB / Hz * np.array(data[1]) / 3.14e13 # kelvin / erg * s * erg / Myr
    plt.plot(data[0],1.4*dT,'g',label=r'$\alpha = 2.5$')

    plt.legend(loc='upper right',fontsize=20)

plt.yscale('log')
plt.xlabel(r'$z$', size=28)
plt.axis()#[1e-3,10,1e-2,1e5],interpolation='none')
plt.xlim([6,20])

#data = read_file('output/test_with_CR_100_keV_igm.txt',0,9)
#plt.plot(data[0],data[1],'r',label='Ph-Heating')

#data = read_file('output/test_with_CR_100_keV_igm.txt',0,10)
#plt.plot(data[0],data[1],'r--',label='CR-Heating')

#plt.text(16.5,10,'Coulomb',size=20,rotation=65)

plot_heating()

#plot_ionization()

#plt.show()

plt.savefig('heating_rate.pdf', format='pdf', dpi=300)