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
    ylabel = r'energy loss length [Mpc]'
    xscale = 'log'
    yscale = 'log'
    xlim = (1, 1e3)
    ylim = (1e-1, 1e4)

    # Configure the axes with the above settings
    set_axes(ax, xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim)

    # Load data
    filename = '../build/diffusion_time_output.txt'
    E, t_H, v, t_a, t_ion, t_C, t_CG = np.loadtxt(filename, usecols=(0, 1, 2, 6, 7, 8, 9), unpack=True)

    # Plot the data
    ax.plot(E, v * t_a, label='adiabatic', color='tab:red')
    ax.plot(E, v * t_ion, label=r'ionization [$x_e = 0$]', color='tab:green')
    ax.plot(E, v * t_C, label='Coulomb [$x_e = 1$]', color='tab:blue')
    ax.plot(E, v * t_CG, label='Coulomb [$x_e = 10^{-3}$]', ls=':', color='tab:blue')

    # Add a legend
    ax.legend(fontsize=20)

    # Save the figure using save_figure function with confirmation message
    save_figure(fig, filepath + 'ccrh_losses_z10.pdf')

if __name__ == "__main__":
    plot_losses()
