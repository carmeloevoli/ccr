import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('ccrh.mplstyle')
import numpy as np

def plot_hmf():
    def set_axes(ax):
        ax.set_xlabel(r'M [M$_\odot$]')
        ax.set_xscale('log')
        ax.set_xlim([1e6, 1e12])
        #ax.set_yscale('log')
        ax.set_ylabel(r'M$^{5/2}$ dN/dM [a.u.]')
        ax.set_ylim([0, 1.1])

    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    filename = '../hmfFrom21cmFast/hmf_allz.txt'

    M, dNdM_10, dNdM_6, dNdM_3, dNdM_1 = np.loadtxt(filename, usecols=(0,1,2,3,4), unpack=True)

    SFR = np.power(M, 1.5)

    y = SFR * M * dNdM_10
    ax.plot(M, y / max(y), color='tab:orange', label='z = 10')

    y = SFR * M * dNdM_6
    ax.plot(M, y / max(y), color='tab:blue', label='z = 6')

    ax.vlines(1.6e9, 0, 1.3, ls=':', color='tab:gray')
    ax.vlines(5.5e10, 0, 1.3, ls=':', color='tab:gray')

    ax.legend()
    plt.savefig('ccrh_hmf_z10.pdf')

def plot_meandistance():
    def set_axes(ax):
        ax.set_xlabel(r'M [M$_\odot$]')
        ax.set_xscale('log')
        ax.set_xlim([1e6, 1e12])
        ax.set_yscale('log')
        ax.set_ylabel(r'd [Mpc]')
        ax.set_ylim([0.1, 1e2])

    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    filename = '../hmfFrom21cmFast/hmf_allz.txt'

    M, dNdM_10, dNdM_6, dNdM_3, dNdM_1 = np.loadtxt(filename, usecols=(0,1,2,3,4), unpack=True)

    y = np.power(M * dNdM_10, -0.34)
    ax.plot(M, y, color='tab:orange', label='z = 10')    

    y = np.power(M * dNdM_6, -0.34)
    ax.plot(M, y, color='tab:blue', label='z = 6')    

    ax.vlines(1.6e9, 0, 1e3, ls=':', color='tab:gray')
    ax.hlines(2.2, 1e6, 1e12, ls=':', color='tab:gray')

    ax.vlines(5.5e10, 0, 1e3, ls=':', color='tab:gray')
    ax.hlines(6.2, 1e6, 1e12, ls=':', color='tab:gray')

    ax.legend()
    plt.savefig('ccrh_meand_z10.pdf')

def plot_timescales():
    def set_axes(ax):
        ax.set_xlabel(r'E [MeV]')
        ax.set_xscale('log')
        ax.set_xlim([0.1, 1e4])
        ax.set_ylabel(r'timescale [Gyr]')
        ax.set_yscale('log')
        ax.set_ylim([1e-3, 1e3])

    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    filename = '../build/timescales_z10.txt'

    E, t_H, t_light, t_D, t_pp, t_a, t_ion, t_C, t_CG = np.loadtxt(filename, usecols=(0,1,2,3,4,5,6,7,8), unpack=True)

    ax.plot(E, t_H, label='Universe age')
    ax.plot(E, t_light, label='ligth')
    ax.plot(E, t_D, label='B = 10$^{-16}$ G')
    ax.plot(E, t_pp, label='pp')
    ax.plot(E, t_a, label='a')
    ax.plot(E, t_ion, label='ion')
    ax.plot(E, t_C, label='C')
    ax.plot(E, t_CG, label='G', ls=':')

    ax.legend()
    plt.savefig('ccrh_timescales_z10.pdf')

if __name__== "__main__":
    #plot_hmf()
    #plot_meandistance()
    plot_timescales()
    #plot_larmor()