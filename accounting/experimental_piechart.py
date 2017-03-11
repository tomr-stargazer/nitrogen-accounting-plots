"""
I'd like to make pie-chart plots for the nitrogen accounting... let's see what we can come up with.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt


def test_piechart():

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    labels = "HCN jk", "CH3CN nope", "something", "yes"
    fractions = [0.45, 0.2, 0.2, 0.15]

    ax.pie(fractions, labels=labels)
    plt.show()

    return fig


