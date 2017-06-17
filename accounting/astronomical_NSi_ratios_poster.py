"""
This is ONLY to make the N/Si ratios for the poster.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

def make_N_Si_figure():

    # Orion_KL_range = np.array([ 0.00510412  ,  0.02398938,  0.1232646 ])
# array([ 0.02398938,  0.1232646 ,  0.00510412])
# 1.51612903e-06,   7.79032258e-06,   3.22580645e-07
# I actually want: 0.0043, mean([0.006, 0.010]), 0.014
    Orion_KL_range = np.array([0.0043, np.mean([0.006, 0.010]), 0.014])
    # IRAS_16293_range = np.array([0.0015, 0.003, 0.023])
    IRAS_16293_range = np.array([0.0042, 0.0058, 0.0082])

#    NGC6334I_range = np.array([0.01387131, 0.00693565, 0.02080696])
# 1.74e-06 2.436e-06 1.044e-06
    NGC6334I_range = np.array([0.02753164556962025, 0.03854430379746835, 0.01651898734177215])    

    fig = plt.figure(figsize=(2, 3), dpi=400)
    ax = fig.add_subplot(111)

    for i, source in enumerate([Orion_KL_range, IRAS_16293_range, NGC6334I_range]):

        min_val = np.min(source)
        max_val = np.max(source)
        median_val = np.median(source)

        print(min_val, median_val, max_val)

        ax.errorbar([i], [median_val], yerr=[[median_val-min_val], [max_val-median_val]], 
            fmt='o', color='tab:red', ms=4, ecolor='tab:red', capsize=4, elinewidth=1, capthick=1)

    ax.yaxis.tick_right()
    ax.tick_params(axis='y', direction='in', which='both', left='on')
    ax.tick_params(axis='x', labelbottom='off', bottom='off')

    plt.semilogy()
    plt.xlim(-0.5,2.5)
    plt.ylim(1e-4, 40)

    # plt.show()
    plt.savefig("OrionKL_NGC_IRAS16293_N_Si_small.png", bbox_inches='tight')
    return fig


