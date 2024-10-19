import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('ccrh.mplstyle')
import numpy as np

from utils import set_axes, save_figure

def plot_spectrum(filepath='figs/'):
    # Create the figure and axis
    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)

    # Set the axes using the generalized set_axes function
    xlabel = r'E [MeV]'
    ylabel = r'E$^2$ Q [10$^{-33}$ erg/cm$^3$/s]'
    xscale = 'log'
    yscale = 'log'
    xlim = (1e0, 1e6)
    ylim = (1e-3, 1)
    
    # Configure the axes with the above settings
    set_axes(ax, xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim)

    # Load data
    filename = '../build/energy_spectrum_output.txt'
    E, Q20, Q24, Q27 = np.loadtxt(filename, usecols=(0, 1, 2, 3), unpack=True)

    # Plot the data
    ax.plot(E, Q20 / 1e-33, label='2.0', color='tab:gray')
    ax.plot(E, Q24 / 1e-33, label='2.4', color='tab:orange')
    ax.plot(E, Q27 / 1e-33, label='2.7', color='tab:blue')

    # Add a legend
    ax.legend(fontsize=20)

    # Save the figure using save_figure function with confirmation message
    save_figure(fig, filepath + 'ccrh_qspectrum.pdf')


def plot_solution(filepath='figs/'):
    # Create the figure and axis
    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)

    # Set the axes using the generalized set_axes function
    xlabel = r'E [MeV]'
    ylabel = r'E$^2$ N'
    xscale = 'log'
    yscale = 'log'
    xlim = (0.1, 1e6)
    ylim = (1e-10, 1e-6)
    
    # Configure the axes with the above settings
    set_axes(ax, xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim)

    # Load data
    filename = '../build/ccrh_solution.txt'
    E, N, Q = np.loadtxt(filename, usecols=(0, 1, 2), unpack=True)
    E2 = E * E

    # Plot the data
    ax.plot(E, E2 * N , label='2.0', color='tab:gray')
    ax.plot(E, E2 * Q , label='2.0', color='tab:gray', ls=':')

    # Save the figure using save_figure function with confirmation message
    save_figure(fig, filepath + 'ccrh_solution.pdf')

def plot_rates(filepath='figs/'):
    # Create the figure and axis
    fig = plt.figure(figsize=(13.5, 8.5))
    ax = fig.add_subplot(111)

    # Set the axes using the generalized set_axes function
    xlabel = r'E$_{\rm min}$ [MeV]'
    ylabel = r'Normalizaed Rates'
    xscale = 'log'
    yscale = 'linear'
    xlim = (0.1, 1e3)
    ylim = (0, 1)
    
    # Configure the axes with the above settings
    set_axes(ax, xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim)

    # Load data
    filename = '../build/ccrh_solution.txt'
    E, I, C = np.loadtxt(filename, usecols=(0, 3, 4), unpack=True)

    ax.plot(E, I / I[0], color='tab:blue', label='ionization')
    ax.fill_between(E, 0., I / I[0], color='tab:blue', alpha=0.1)
    
    ax.plot(E, C / C[0], color='tab:orange', label='Coulomb')
    ax.fill_between(E, 0., C / C[0], color='tab:orange', alpha=0.1)

    ax.hlines(0.1, 0.1, 1e3, ls=':', color='tab:gray')

    # Save the figure using save_figure function with confirmation message
    save_figure(fig, filepath + 'ccrh_rates.pdf')

if __name__ == "__main__":
    plot_spectrum()
    plot_solution()
    plot_rates()