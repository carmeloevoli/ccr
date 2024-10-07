import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('ccrh.mplstyle')
import numpy as np

from utils import set_axes, save_figure

def plot_timescales(filepath='figs/'):
    # Create the figure and axis
    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)

    # Set the axes using the generalized set_axes function
    xlabel = r'E [MeV]'
    ylabel = r'timescale [Gyr]'
    xscale = 'log'
    yscale = 'log'
    xlim = (1, 1e4)
    ylim = (1e-4, 1e1)
    
    # Configure the axes with the above settings
    set_axes(ax, xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim)

    # Load data
    filename = '../build/diffusion_time_output.txt'
    E, t_H, t_light, t_D = np.loadtxt(filename, usecols=(0, 1, 3, 4), unpack=True)

    # Plot the data
    z = 10.
    ax.plot(E, t_H, label='Universe age at z = 10', color='tab:gray')
    ax.plot(E, t_light, label=r'$\langle d \rangle$ / v', color='tab:orange')
    ax.fill_between(E, 1e-5, t_light, color='tab:orange', alpha=0.2)
    ax.plot(E, t_D, label='B = 10$^{-16}$ G', color='tab:red')

    # Add a legend
    ax.legend(fontsize=20)

    # Save the figure using save_figure function with confirmation message
    save_figure(fig, filepath + 'ccrh_timescales_z10.pdf')

if __name__ == "__main__":
    plot_timescales()
