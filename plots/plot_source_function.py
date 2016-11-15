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
rcParams['lines.linewidth'] = 5
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

def read_file(datafile,xcol,ycol,zcol):
    x = []
    y = []
    z = []
    f = open(datafile,'r')
    header = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        x.append(float(columns[xcol]))
        y.append(float(columns[ycol]))
        z.append(float(columns[zcol]))
    f.close()
    data=[]
    data.append(np.array(x))
    data.append(np.array(y))
    data.append(np.array(z))
    return data

plt.yscale('log')
plt.xscale('log')

plt.title(r'$z=10$',fontsize=32)
plt.xlabel(r'$E$ [GeV]',fontsize=32)
plt.ylabel(r'$E^2$ q$_{\rm p}$ [erg cm$^{-3}$ s$^{-1}$]',fontsize=32)

plt.axis()#[1e7,1e11,1e3,1e9],interpolation='none')

GeV = 1. / 624.151

data = read_file('output/test_new_2_10_spectrum.txt',0,1,2)

plt.plot(data[0] / GeV, data[0]**2.0 * data[1] * data[2], 'b', label=r'$\alpha = 2$')

data = read_file('output/test_new_2.2_10_spectrum.txt',0,1,2)

plt.plot(data[0] / GeV, data[0]**2.0 * data[1] * data[2], 'g', label=r'$\alpha = 2.2$')

data = read_file('output/test_new_2.5_10_spectrum.txt',0,1,2)

plt.plot(data[0] / GeV, data[0]**2.0 * data[1] * data[2], 'r', label=r'$\alpha = 2.5$')

plt.legend(loc='lower right',frameon=False,fontsize=26)

plt.xlim([1e-3,1e2])
plt.ylim([5e-36,.2e-32])

plt.savefig('source_function_z10.pdf',format='pdf',dpi=300)

#plt.show()
