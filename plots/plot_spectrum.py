import matplotlib.pyplot as plt
import numpy as np

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

plt.xlabel(r'$E$ [GeV]',fontsize=16)
plt.ylabel(r'$E^2 CR \, []$',fontsize=16)

plt.axis()#[1e7,1e11,1e3,1e9],interpolation='none')

#plt.xlim([0.1,1])

alpha = 1.

data = read_file('output/test_spectra_192_at_29_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0],data[0]**alpha*data[1],'b')

data = read_file('output/test_spectra_192_at_20_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0],data[0]**alpha*data[1],'r')

data = read_file('output/test_spectra_192_at_10_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0],data[0]**alpha*data[1],'g')

data = read_file('output/test_spectra_192_at_1_fesc_0.0002_fsfr_0.04.txt',0,1)
plt.plot(data[0],data[0]**alpha*data[1],'m')

data = read_file('output/spectra_no_192_at_29.txt',0,1)
#plt.plot(data[0],data[0]**alpha*data[1],'b:')

data = read_file('output/spectra_no_192_at_28.txt',0,1)
#plt.plot(data[0],data[0]**alpha*data[1],'r:')

data = read_file('output/spectra_no_192_at_27.txt',0,1)
#plt.plot(data[0],data[0]**alpha*data[1],'g:')

#plt.legend(loc='upper right')

#plt.savefig('hmf_SFR.pdf',format='pdf',dpi=300)

plt.show()