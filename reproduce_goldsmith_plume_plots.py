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


def make_goldsmith_figure_10(save=True):
    """ Generates my version of Fig. 10 from Goldsmith+'97 """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_array = np.logspace(1, 7, 100)
    temperature_list = [10, 30, 50]

    f_c_array_10K = np.zeros_like(density_array)
    f_c_array_30K = np.zeros_like(density_array)
    f_c_array_50K = np.zeros_like(density_array)

    f_c_array_list = [f_c_array_10K, f_c_array_30K, f_c_array_50K]

    for j, temp in enumerate(temperature_list):

        # Here we have to specify 2 of 3: abundance, temperature, density.
        R = pyradex.Radex(species='c18o', abundance=1.7e-7, column=1e14, temperature=temp)

        for i, rho in enumerate(density_array):

            new_T = R(collider_densities={'H2':rho})

            # use the correction_factor code here!
            f_c = correction_factor_given_radex_table(new_T, [0, 1])

            f_c_array_list[j][i] = f_c
            
    ax.plot(np.log10(density_array), f_c_array_10K, '-', lw=1.5)
    ax.plot(np.log10(density_array), f_c_array_30K, ':', lw=2)
    ax.plot(np.log10(density_array), f_c_array_50K, '--', lw=1.5)

    ax.set_ylim(1, 3)
    ax.set_xlim(2, 7.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("CORRECTION FACTOR FOR COLUMN DENSITY", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(4.5, 2.8, "C$^{18}$O", fontsize=18, family='serif')
    ax.text(6, 2.8, r"$T_k = 50\ \rm{K}$", fontsize=16, family='serif', color='r')
    ax.text(6, 2.2, r"$T_k = 30\ \rm{K}$", fontsize=16, family='serif', color='g')
    ax.text(6, 1.55, r"$T_k = 10\ \rm{K}$", fontsize=16, family='serif', color='b')
    ax.text(2.5, 2.5, "$J=1-0$ & $J=2-1$ OBSERVED", fontsize=14, family='serif')

    ax.minorticks_on()

    if save:
        fig.savefig("goldsmith97_fig_10_reproduction.pdf")

    plt.show()

    return fig


def make_goldsmith_figure_11(save=True):
    """ Generates my version of Fig. 11 from Goldsmith+'97 """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_array = np.logspace(1, 7, 100)
    temperature_list = [10, 30, 50]

    f_c_array_10K = np.zeros_like(density_array)
    f_c_array_30K = np.zeros_like(density_array)
    f_c_array_50K = np.zeros_like(density_array)

    f_c_array_list = [f_c_array_10K, f_c_array_30K, f_c_array_50K]

    for j, temp in enumerate(temperature_list):

        # Here we have to specify 2 of 3: abundance, temperature, density.
        R = pyradex.Radex(species='c18o', abundance=1.7e-7, column=1e14, temperature=temp)

        for i, rho in enumerate(density_array):

            new_T = R(collider_densities={'H2':rho})

            # use the correction_factor code here!
            f_c = correction_factor_given_radex_table(new_T, [0, 1, 2])

            f_c_array_list[j][i] = f_c

            
    ax.plot(np.log10(density_array), f_c_array_10K, '-', lw=1.5)
    ax.plot(np.log10(density_array), f_c_array_30K, ':', lw=2)
    ax.plot(np.log10(density_array), f_c_array_50K, '--', lw=1.5)

    ax.set_ylim(1, 3)
    ax.set_xlim(2, 7.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("CORRECTION FACTOR FOR COLUMN DENSITY", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(4.5, 2.8, "C$^{18}$O", fontsize=18, family='serif')
    ax.text(6, 2, r"$T_k = 50\ \rm{K}$", fontsize=16, family='serif', color='r')
    ax.text(6, 1.55, r"$T_k = 30\ \rm{K}$", fontsize=16, family='serif', color='g')
    ax.text(6, 1.2, r"$T_k = 10\ \rm{K}$", fontsize=16, family='serif', color='b')
    ax.text(2.5, 2.5, "$J=1-0$, $J=2-1$, & $J=3-2$ OBSERVED", fontsize=14, family='serif')

    ax.minorticks_on()

    if save:
        fig.savefig("goldsmith97_fig_11_reproduction.pdf")

    plt.show()

    return fig