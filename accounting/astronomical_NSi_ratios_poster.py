"""
This is ONLY to make the N/Si ratios for the poster.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

def make_N_Si_figure():

    Orion_KL_range = np.array([ 0.009902  ,  0.02398938,  0.1232646 ])

    IRAS_16293_range = np.array([0.0015, 0.003, 0.023])

    NGC6334I_range = np.array([0.01387131, 0.00693565, 0.02080696])

    fig = plt.figure(figsize=(1.5, 3), dpi=400)
    ax = fig.add_subplot(111)

    for i, source in enumerate([Orion_KL_range, IRAS_16293_range, NGC6334I_range]):

        min_val = np.min(source)
        max_val = np.max(source)
        median_val = np.median(source)

        print min_val, median_val, max_val

        ax.errorbar([i], [median_val], yerr=[[median_val-min_val], [max_val-median_val]], 
            fmt='o', color='tab:red', ecolor='tab:red', capsize=4, elinewidth=2, capthick=2)

    plt.semilogy()
    plt.xlim(-1,3)
    plt.ylim(1e-4, 1)

    plt.show()
    plt.savefig("OrionKL_NGC_IRAS16293_N_Si_small.png", bbox_inches='tight')
    return fig


