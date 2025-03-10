import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.cm as cm

def plot_series(dependent_vars, independent_var, filename=None):
    """
    Plots multiple dependent variables against a single independent variable.
    
    Parameters:
    dependent_vars (dict): Dictionary mapping series names to dependent variable arrays.
    independent_var (dict): Dictionary of length 1 mapping a series name to an independent variable array.
    filename (str, optional): If provided, saves the plot as an SVG file with the given filename.
    """
    if len(independent_var) != 1:
        raise ValueError("independent_var must contain exactly one key-value pair.")
    
    x_label, x_values = list(independent_var.items())[0]  # Extract the single independent variable
    
    num_lines = len(dependent_vars)
    inferno_colors = cm.inferno(np.linspace(0, 1, num_lines))  # Generate colors from inferno colormap
    
    fig, ax = plt.subplots()
    
    handles = []
    for (label, y_values), color in zip(dependent_vars.items(), inferno_colors):
        line, = ax.plot(x_values, y_values, label=label, color=color)  # No markers
        handles.append(line)
    
    # Create a discrete legend
    ax.legend(handles=handles, labels=dependent_vars.keys(), loc='best', frameon=True)
    
    ax.set_xlabel(x_label)  # Label x-axis
    ax.set_ylabel("")  # No title on y-axis
    
    ax.yaxis.set_ticks_position('left')  # Enable y-axis ticks
    ax.xaxis.set_ticks_position('bottom')
    
    ax.grid(True)  # Enable gridmarks
    
    if filename:
        plt.savefig(filename, format='svg')  # Save as SVG if filename is provided
    
    plt.show()
