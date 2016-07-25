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

plt.xticks(size=30)
plt.yticks(size=30)
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
plt.xscale('log')

plt.xlabel(r'$E$ [GeV]',fontsize=30)
plt.ylabel(r'$E^2$ N$_p$ [erg]',fontsize=30)

plt.axis()#[1e7,1e11,1e3,1e9],interpolation='none')

alpha = 2.
GeV = 1. / 624.151

data = read_file('output/test_with_CR_keV_spectra_224_at_20_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0]/GeV,data[0]**alpha*data[1],'r',label=r'$z=20$')

data = read_file('output/test_with_CR_no_losses_keV_spectra_224_at_20_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0]/GeV,data[0]**alpha*data[1],'r--')

data = read_file('output/test_with_CR_keV_spectra_224_at_10_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0]/GeV,data[0]**alpha*data[1],'g',label=r'$z=10$')

data = read_file('output/test_with_CR_no_losses_keV_spectra_224_at_10_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0]/GeV,data[0]**alpha*data[1],'g--')

data = read_file('output/test_with_CR_keV_spectra_224_at_6_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0]/GeV,data[0]**alpha*data[1],'b',label=r'$z=6$')

data = read_file('output/test_with_CR_no_losses_keV_spectra_224_at_6_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0]/GeV,data[0]**alpha*data[1],'b--')

plt.legend(loc='upper right',fontsize=24)

plt.xlim([1e-4,1e0])
plt.ylim([1e-20,1e-12])

plt.savefig('proton_evolution.pdf',format='pdf',dpi=300)

#plt.show()