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

def read_SFR(datafile):
    x0 = []
    x1 = []
    x2 = []
    x3 = []
    x4 = []
    f = open(datafile,'r')
    header = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        x0.append(float(columns[0]))
        x1.append(float(columns[1]))
        x2.append(float(columns[2]))
        x3.append(float(columns[3]))
        x4.append(float(columns[4]))
    f.close()
    data=[]
    data.append(np.array(x0))
    data.append(np.array(x1))
    data.append(np.array(x2))
    data.append(np.array(x3))
    data.append(np.array(x4))
    return data

def plot_SFR(data,color):
    z = 0.5 * (data[1] + data[0])
    dz_lo = z - data[0]
    sfr = data[2]
    sfr_lo = data[3]
    plt.errorbar(z, sfr, xerr=dz_lo, yerr=sfr_lo, fmt='D', markersize='8', elinewidth=2, capsize=6, capthick=2, color=color, label='Madau \& Dickinson 2014')

def plot_SFR_theory(filename,fs,color,label):
    mass_sun = 1.989e33
    Mpc = 3.0856775807e24
    yr = 3.14e7
    data = read_output_file(filename,0,2)
    plt.plot(data[0],np.log10(fs * data[1] / (mass_sun / Mpc**3 / yr)), color = color)


#data = read_SFR("SFR_IR_Madau14.txt")

#plot_SFR(data,'k')

data = read_SFR("SFR_UV_Madau14.txt")

plot_SFR(data,'k')

fs = 0.01

plot_SFR_theory('ccr/SFR_new_60.txt',fs,'g','$V_c = 30$ km/s')

fs = 0.012

plot_SFR_theory('ccr/SFR_new_100.txt',fs,'b','$V_c = 50$ km/s')

fs = 0.021

plot_SFR_theory('ccr/SFR_new_200.txt',fs,'r','$V_c = 100$ km/s')

#plot_SFR_theory('ccr/SFR.txt',fs,'m','$V_c = 0$ km/s')

plt.xscale('log')

x = [1,2,3,4,5,6,7,8,9,10]

labels = ['1','2','3','4','5','6','7','8','9','10']

plt.xticks(x, labels)

plt.xlabel(r'z', size=30)
plt.ylabel(r'log SFR [M$_\odot$ Mpc$^{-3}$ yr$^{-1}$]', size=30)

plt.xlim([1,10])
plt.ylim([-2.5,-0.5])

plt.legend(fontsize=22,loc='lower left',frameon=False)

plt.savefig('SFR_vs_Madau.pdf', format='pdf', bbox_inches='tight', dpi=300)

#plt.show()
