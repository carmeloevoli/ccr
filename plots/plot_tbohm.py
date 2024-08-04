import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('ccrh.mplstyle')
import numpy as np

def plot_hmf():
    def min_mass(z):
        return 1e8 * (10. / (1. + z))**(1.5)

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

    m = min_mass(10.)
    ax.fill_between([1e6, m], 0, 1.1, color='tab:blue', alpha=0.3, label='No SFR')

#    y = SFR * M * dNdM_6
#    ax.plot(M, y / max(y), color='tab:blue', label='z = 6')

    ax.vlines(1.6e9, 0, 1.3, ls=':', color='tab:gray')
#    ax.vlines(5.5e10, 0, 1.3, ls=':', color='tab:gray')

    ax.legend()
    plt.savefig('ccrh_hmf_z10.pdf')

def plot_meandistance():
    def set_axes(ax):
        ax.set_xlabel(r'M [M$_\odot$]')
        ax.set_xscale('log')
        ax.set_xlim([1e7, 1e12])
        ax.set_yscale('log')
        ax.set_ylabel(r'$\langle d \rangle$ [Mpc]')
        ax.set_ylim([0.1, 1e2])

    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    filename = '../hmfFrom21cmFast/hmf_allz.txt'

    M, dNdM_10, dNdM_6, dNdM_3, dNdM_1 = np.loadtxt(filename, usecols=(0,1,2,3,4), unpack=True)

    y = np.power(M * dNdM_10, -0.34)
    ax.plot(M, y, color='tab:orange', label='z = 10')    

#    y = np.power(M * dNdM_6, -0.34)
#    ax.plot(M, y, color='tab:blue', label='z = 6')    

    ax.vlines(1.6e9, 0, 1e3, ls=':', color='tab:gray')
    ax.hlines(2.2, 1e6, 1e12, ls=':', color='tab:gray')

#    ax.vlines(5.5e10, 0, 1e3, ls=':', color='tab:gray')
#    ax.hlines(6.2, 1e6, 1e12, ls=':', color='tab:gray')

    ax.legend()
    plt.savefig('ccrh_meand_z10.pdf')

def plot_timescales():
    def set_axes(ax):
        ax.set_xlabel(r'E [MeV]')
        ax.set_xscale('log')
        ax.set_xlim([1, 1e4])
        ax.set_ylabel(r'timescale [Gyr]')
        ax.set_yscale('log')
        ax.set_ylim([1e-3, 1e2])

    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    filename = '../build/timescales_z10.txt'

    E, t_H, t_light, t_D, t_pp, t_a, t_ion, t_C, t_CG = np.loadtxt(filename, usecols=(0,1,2,3,4,5,6,7,8), unpack=True)

    ax.plot(E, t_H, label='Universe age at z = 10', color='tab:gray')
    ax.plot(E, t_light, label=r'$\langle d \rangle$ / v', color='tab:orange')
    ax.fill_between(E, 1e-3, t_light, color='tab:orange', alpha=0.2)
    ax.plot(E, t_D, label='B = 10$^{-16}$ G', color='tab:red')
    
    ax.legend(fontsize=20)
    plt.savefig('ccrh_timescales_z10.pdf')

def plot_losses():
    def set_axes(ax):
        ax.set_xlabel(r'E [MeV]')
        ax.set_xscale('log')
        ax.set_xlim([1, 1e4])
        ax.set_ylabel(r'timescale [Gyr]')
        ax.set_yscale('log')
        ax.set_ylim([1e-3, 1e2])
        
    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    filename = '../build/timescales_z10.txt'

    E, t_H, t_light, t_D, t_pp, t_a, t_ion, t_C, t_CG = np.loadtxt(filename, usecols=(0,1,2,3,4,5,6,7,8), unpack=True)

    ax.plot(E, t_H, label='Universe age at z = 10', color='tab:gray')
    ax.plot(E, t_light, label=r'$\langle d \rangle$ / v', color='tab:orange')
    ax.fill_between(E, 1e-3, t_light, color='tab:orange', alpha=0.2)
    
    ax.plot(E, t_a, label='adiabatic', color='tab:red')
    ax.plot(E, t_ion, label=r'ionization [$x_e = 0$]', color='tab:green')
    ax.plot(E, t_C, label='Coulomb [$x_e = 1$]', color='tab:blue')
    ax.plot(E, t_CG, label='Coulomb [$x_e = 10^{-3}$]', ls=':', color='tab:blue')

    #ax.plot(E, t_D, label='B = 10$^{-16}$ G', color='tab:red')
    #    ax.plot(E, t_pp, label='pp')

    ax.legend(fontsize=20)
    plt.savefig('ccrh_losses_z10.pdf')

def plot_spectrum():
    def spectrum(E, slope):
        mpc2 = 0.938e3
        p = np.sqrt(E * E + 2. * mpc2 * E)
        beta = p / (E + mpc2)
        p_0 = 1e3
        return np.power(p / p_0, 3. - slope)
        
    def set_axes(ax):
        ax.set_xlabel(r'E [MeV]')
        ax.set_xscale('log')
        ax.set_xlim([1, 1e4])
        ax.set_ylabel(r'$E_{\rm CR}$')
        #ax.set_yscale('log')
        ax.set_ylim([1e2, 5e2])

    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    E = np.logspace(0, 10, 1000)
    
    f = spectrum(E, 4.5)
    
    ax.plot(E, E * f)
    
    plt.savefig('ccrh_spectrum.pdf')

#// double spectrum(const double& E_k, const double& SN_slope) {
#//     double p = sqrt(E_k * E_k + 2. * mass_proton_c2 * E_k);
#//     double beta = p / (E_k + mass_proton_c2);
#//     double E_k_0 = reference_energy;
#//     double p_0 = sqrt(E_k_0 * E_k_0 + 2. * mass_proton_c2 * E_k_0);
#
#//     return 1. / beta * pow(p / p_0, -SN_slope);
#// }


if __name__== "__main__":
#    plot_hmf()
#    plot_meandistance()
#    plot_timescales()
#    plot_losses()
    plot_spectrum()
    #TBD Larmor!
