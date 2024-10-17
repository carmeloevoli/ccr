import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('ccrh.mplstyle')
import numpy as np

from utils import set_axes, save_figure

def plot_losses(filepath='figs/'):
    # Create the figure and axis
    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)

    # Set the axes using the generalized set_axes function
    xlabel = r'E [MeV]'
    ylabel = r'energy loss relative timescale'
    xscale = 'log'
    yscale = 'log'
    xlim = (1e-2, 1e4)
    ylim = (1e-2, 1e2)

    # Configure the axes with the above settings
    set_axes(ax, xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim)

    # Load data
    filename = '../build/losses_timescales_output.txt'
    E, t_C_in, t_C_out, t_ion, t_a = np.loadtxt(filename, usecols=(0, 1, 2, 3, 4), unpack=True)

    # Plot the data
    ax.plot(E, t_ion / t_a, label=r'ionization [$x_e = 0$]', color='tab:blue')
    ax.plot(E, t_C_out / t_a, label='Coulomb [$x_e = 1$]', color='tab:orange')
    ax.plot(E, t_C_in / t_a, label='Coulomb [$x_e = 10^{-3}$]', ls='-.', color='tab:orange')

    ax.hlines(1., 1e-2, 1e8, ls=':', color='tab:gray')

    ax.vlines(1.08112e-1, 1e-3, 1e2, ls='-', color='tab:orange')
    ax.vlines(5.73160e+0, 1e-3, 1e2, ls='-', color='tab:blue')
    ax.vlines(8.77068e+0, 1e-3, 1e2, ls='-.', color='tab:orange')

    ax.text(1e3, 12, 'z = 12')
    # Add a legend
    ax.legend(fontsize=20)

    # Save the figure using save_figure function with confirmation message
    save_figure(fig, filepath + 'ccrh_losses_z10.pdf')

if __name__ == "__main__":
    plot_losses()
