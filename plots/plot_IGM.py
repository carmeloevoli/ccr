#!/bin/bash/python
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import numpy as np

# BEGIN PREAMBLE

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Helvetica Neue')
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)
rcParams['legend.numpoints'] = 1
rcParams['lines.linewidth'] = 4
rcParams['figure.autolayout'] = True

fig = plt.figure(figsize=(8.3, 7.8))
ax = fig.add_subplot(1, 1, 1)

for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)

ax.minorticks_on()
ax.tick_params('both', length=15, width=1.5, which='major', pad=6)
ax.tick_params('both', length=10, width=1.3, which='minor', pad=6)

plt.xticks(size=30)
plt.yticks(size=30)

# END PREAMBLE

def read_output_file(filename,xcol,ycol):
    x0 = []
    x1 = []
    f = open(filename,'r')
    header = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        x0.append(float(columns[xcol]))
        x1.append(float(columns[ycol]))
    f.close()
    data=[]
    data.append(np.array(x0))
    data.append(np.array(x1))
    return data

data = read_output_file('output/test_new_igm.txt',0,1)

ax1 = plt.subplot(211)
ax1.minorticks_on()
ax1.tick_params('both', length=15, width=1.5, which='major', pad=6)
ax1.tick_params('both', length=10, width=1.3, which='minor', pad=6)

plt.xlim([0,20])
plt.ylim([1e-4,1e0])
plt.yscale('log')
plt.plot(data[0],data[1],label=r'$x_{\rm HII}$')
plt.plot(data[0],1.0-data[1],label=r'$x_{\rm HI}$')
plt.setp(ax1.get_yticklabels(), fontsize=26)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.legend(loc=1,frameon=False,fontsize=20)

plt.ylabel(r'$x$', size=26)

data = read_output_file('output/test_new_igm.txt',0,3)

ax2 = plt.subplot(212, sharex=ax1)
#ax2.minorticks_on()
ax2.tick_params('both', length=15, width=1.5, which='major', pad=6)
ax2.tick_params('both', length=10, width=1.3, which='minor', pad=6)

plt.xlim([0,20])
plt.ylim([0,0.07])
plt.plot(data[0],data[1])
plt.setp(ax2.get_yticklabels(), fontsize=26)
plt.setp(ax2.get_xticklabels(), fontsize=26)

plt.errorbar([0.1,0.1], [0.055, 0.055], yerr=[0.009, 0.009], fmt='D', markersize='8', elinewidth=2, capsize=6, capthick=2, color='r')

plt.text(1,0.06,'PLANCK',fontsize=20,color='r')
plt.text(13,0.05,r'$f_{\rm esc} = 10^{-2}$',fontsize=25)

plt.xlabel(r'z', size=26)
plt.ylabel(r'$\tau_e$', size=26)

plt.savefig('IGM_evolution.pdf', format='pdf', bbox_inches='tight', dpi=300)
