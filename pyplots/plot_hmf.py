import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('ccrh.mplstyle')
import numpy as np

from utils import set_axes, save_figure

# Constants
LINE_COLOR = 'tab:gray'
SFR_COLOR = 'tab:orange'
NO_SFR_COLOR = 'tab:blue'

def min_mass(z):
    """Calculate the minimum mass as a function of redshift."""
    return 1e8 * (10. / (1. + z)) ** 1.5

def load_data(filename):
    """Load mass function data from file."""
    return np.loadtxt(filename, usecols=(0, 1, 2, 3, 4), unpack=True)

def plot_sfr(ax, M, dNdM, z, SFR):
    """Plot star formation rate (SFR) with respect to mass."""
    y = SFR * M * dNdM
    ax.plot(M, y / max(y), color=SFR_COLOR, label=f'z = {z}')

def fill_no_sfr_region(ax, m_min, xlim, ylim):
    """Fill the 'No SFR' region in the plot."""
    ax.fill_between([xlim[0], m_min], ylim[0], ylim[1], 
                    color=NO_SFR_COLOR, alpha=0.3, label='No SFR')

def add_vlines(ax, positions, limits):
    """Add vertical lines at specified positions."""
    for pos in positions:
        ax.vlines(pos, limits[0], limits[1], ls=':', color=LINE_COLOR)

def plot_hmf(filepath='figs/'):
    """Main function to plot the halo mass function."""
    # Set up figure and axis
    fig, ax = plt.subplots(figsize=(13.5, 8.5))
    
    # Generalized axis configuration
    xlabel = r'M [M$_\odot$]'
    ylabel = r'M$^{5/2}$ dN/dM [a.u.]'
    xscale = 'log'
    xlim = (1e6, 1e12)
    ylim = (0, 1.1)
    
    # Set the axes with optional parameters
    set_axes(ax, xlabel=xlabel, ylabel=ylabel, xscale=xscale, xlim=xlim, ylim=ylim)

    # Load data
    filename = '../hmfFrom21cmFast/hmf_allz.txt'
    M, dNdM_10, dNdM_6, dNdM_3, dNdM_1 = load_data(filename)

    # Calculate Star Formation Rate (SFR)
    SFR = np.power(M, 1.5)

    # Plot for z = 10
    plot_sfr(ax, M, dNdM_10, 10, SFR)
    fill_no_sfr_region(ax, min_mass(10), xlim, ylim)

    # Add vertical lines at specific mass values
    add_vlines(ax, [1.6e9], [0, 1.3])

    # Add legend
    ax.legend()

    # Save the figure and print confirmation
    save_figure(fig, filepath + 'ccrh_hmf_z10.pdf')

if __name__ == "__main__":
    plot_hmf()
