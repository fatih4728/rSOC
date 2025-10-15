# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 12:58:52 2025

Here I want to create my ideal 2D-XY-Graphic

@author: smfadurm
"""


# %%
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def standardized_plot(
    x, y_list,
    xlabel="X", ylabel="Y", title=None,
    grid=False, style="seaborn-v0_8-paper",
    savepath=None, labels=None, startAtZero=True,
    font="Helvetica"
):
    """
    x : list or array
    y_list : list of lists/arrays (each inner list/array is a separate line to plot)
    labels : list of strings for legend (optional)
    startAtZero : if True, axes start at 0 with a small offset to avoid tick artifacts
    """
    # Ensure y_list is a list of 1D arrays
    y_list_corrected = []
    for y in y_list:
        y = np.array(y)
        if y.ndim == 2 and y.shape[0] == 1:
            y = y.flatten()
        elif y.ndim > 1 and y.shape[0] != x.shape[0]:
            y = y.T if y.shape[1] == x.shape[0] else y.flatten()
        y_list_corrected.append(y)

    # Use a consistent Matplotlib style
    plt.style.use('default')
    plt.style.use(style)
    
    if font=='Helvetica':
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.sans-serif'] = [font]
    else:
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = [font]

    # Create a new figure and axis
    fig, ax = plt.subplots(figsize=(6*1.3, 4*1.3))
  
    # Plot each dataset and keep handles for legend
    color=['#F46A25', '#87CEEB', 'purple', '#B2B2B2', 'grey', 'lime']

    
    lines = []
    for idx, y in enumerate(y_list_corrected):
        label = labels[idx] if labels and idx < len(labels) else f"Line {idx+1}"
        line, = ax.plot(x, y, linewidth=2, label=label, color=color[idx])
        lines.append(line)
    
    # Label axes and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    
    # Optional grid
    if grid:
        ax.grid(True, linestyle='--', alpha=0.7)
    
    # Configure major and minor ticks: inward-facing on all sides
    ax.tick_params(which='major', direction='in', 
                   top=True, right=True, length=6, width=0.8)
    ax.tick_params(which='minor', direction='in', 
                   top=True, right=True, length=4, width=0.8)

    # Automatically adjust tick intervals with alternating major and minor ticks
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    # Set axes starting at 0 with a small offset to avoid tick artifacts
    if startAtZero:
        xmin = min(x)
        y_max = max([y.max() for y in y_list_corrected])
        y_min = min([y.min() for y in y_list_corrected])
        
        ax.set_xlim(left=xmin, right=max(x))
        ax.set_ylim(bottom=y_min, top=y_max)

    # Add legend
    ax.legend()
    
    # Tidy layout
    fig.tight_layout()
    
    # Save as SVG (vector graphic)
    if savepath:
        if not savepath.lower().endswith(".svg"):
            savepath += ".svg"
        fig.savefig(savepath, format="svg")
    
    return fig, ax




# %%
if __name__=="__main__":
    
    # data
    N = 100
    x = np.linspace(0, 10, N)
    y = np.linspace(1, 100, N)
    y2 = 0.009*y**2
    y3 = y*0.5 + 30
    y4 = y**0.5
    y5 = y4+40
    
    # function call
    standardized_plot(x, [y, y2, y3, y4], 
                      'Distance / m', 'Height / cm',   
                      labels=["base", "power of two", "linear", "sub"],
                      style='seaborn-v0_8-poster',
                      savepath=r'C:\Users\smfadurm\Desktop\test',
                      font = 'Helvetica')














