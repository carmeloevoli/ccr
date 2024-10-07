def set_axes(ax, xlabel=None, ylabel=None, xscale='linear', yscale='linear', xlim=None, ylim=None):
    """
    Configure the axes of the plot with optional parameters.
    
    Parameters:
    - ax: The axis to configure.
    - xlabel: Label for the x-axis (optional).
    - ylabel: Label for the y-axis (optional).
    - xscale: Scale for the x-axis ('linear', 'log', etc., default is 'linear').
    - yscale: Scale for the y-axis ('linear', 'log', etc., default is 'linear').
    - xlim: Tuple specifying x-axis limits (optional).
    - ylim: Tuple specifying y-axis limits (optional).
    """
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

def save_figure(fig, filename):
    """
    Save the figure to the specified path and print a confirmation message.
    
    Parameters:
    - fig: The figure object to save.
    - filename: The figure filename.
    """
    fig.savefig(filename)
    print(f"Figure saved as {filename}")