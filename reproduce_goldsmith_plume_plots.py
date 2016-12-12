"""
This is a 'sanity check' to ensure that I know how to use pyradex 
and how to compute f_c correctly (under benchmark-like conditions).

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import pyradex

# This import will have to change if I refactor correction_factor much.
from correction_factor.correction_factor_2 import correction_factor_given_radex_table

# Default values come from Goldsmith/Bergin/Lis '97, esp. the c18o abundance
def make_goldsmith_figure_9(save=True, column=1e16, abundance=1.7e-7):
    """ Generates my version of Fig. 9 from Goldsmith+'97 """

    # Here we have to specify 2 of 3: abundance, temperature, density.
    R = pyradex.Radex(species='c18o', abundance=abundance, temperature=20, column=column)

    fig = plt.figure()

    ax = fig.add_subplot(111)

    density_array = np.logspace(1, 6.5, 100)

    Tex_array_1_0 = np.zeros_like(density_array)
    Tex_array_2_1 = np.zeros_like(density_array)
    Tex_array_3_2 = np.zeros_like(density_array)

    for i, rho in enumerate(density_array):

        new_T = R(collider_densities={'H2':rho})

        Tex_array_1_0[i] = new_T[0]['Tex']
        Tex_array_2_1[i] = new_T[1]['Tex']
        Tex_array_3_2[i] = new_T[2]['Tex']

    ax.plot(np.log10(density_array), Tex_array_1_0, '-', lw=1.5)
    ax.plot(np.log10(density_array), Tex_array_2_1, ':', lw=2)
    ax.plot(np.log10(density_array), Tex_array_3_2, '--', lw=1.5)

    ax.set_ylim(0,25)
    ax.set_xlim(1,6.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("EXCITATION TEMPERATURE/K", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(1.5, 20, "C$^{18}$O", fontsize=18, family='serif')
    ax.text(2.35, 15, "$J = 1 - 0$", fontsize=16, family='serif', color='b')
    ax.text(3.6, 18, "$J = 2 - 1$", fontsize=16, family='serif', color='g')
    ax.text(4.1, 10, "$J = 3 - 2$", fontsize=16, family='serif', color='r')

    ax.minorticks_on()

    if save:
        fig.savefig("goldsmith97_fig_9_reproduction.pdf")

    plt.show()

    return fig

