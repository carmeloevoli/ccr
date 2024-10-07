import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('ccrh.mplstyle')
import numpy as np

from utils import set_axes, save_figure

def plot_meandistance(filepath='figs/'):
    # Create the figure and axis
    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)

    # Set the axes using the generalized set_axes function
    xlabel = r'M [M$_\odot$]'
    ylabel = r'Mean Distance [a.u.]'
    xscale = 'log'
    yscale = 'log'
    xlim = (1e7, 1e12)
    ylim = (0.1, 1e2)
    
    # Configure the axes with the above settings
    set_axes(ax, xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim)

    # Load data
    filename = '../hmfFrom21cmFast/hmf_allz.txt'
    M, dNdM_10, dNdM_6, dNdM_3, dNdM_1 = np.loadtxt(filename, usecols=(0, 1, 2, 3, 4), unpack=True)

    # Plot the data for z = 10
    y = np.power(M * dNdM_10, -0.34)
    ax.plot(M, y, color='tab:orange', label='z = 10')

    # Add vertical and horizontal lines
    ax.vlines(1.6e9, 0, 1e3, ls=':', color='tab:gray')
    ax.hlines(2.2, 1e6, 1e12, ls=':', color='tab:gray')

    # Add a legend
    ax.legend()

    # Save the figure using save_figure function with confirmation message
    save_figure(fig, filepath + 'ccrh_meand_z10.pdf')

if __name__ == "__main__":
    plot_meandistance()
