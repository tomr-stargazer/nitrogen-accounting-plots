"""
This is my first attempt at making an exploratory N/Si plot.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import ScalarFormatter

def make_fig():

    fig = plt.figure(figsize=(8,6), dpi=100)
    ax = fig.add_subplot(111)

    ax.plot()
    ax.set_xlabel("Arbitrary X axis")
    ax.set_ylabel("N abundance relative to Si", fontsize=16)
    ax.semilogy()
    ax.yaxis.set_major_formatter(ScalarFormatter())

    ax.plot([0.01,100], [0.01, 100], 'ro')

    return fig