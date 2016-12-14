"""
This is a 'sanity check' to ensure that I know how to use pyradex 
and how to compute f_c correctly (under benchmark-like conditions).

"""

from __future__ import division
import time

import numpy as np
import matplotlib.pyplot as plt
import pyradex

# This import will have to change if I refactor correction_factor much.
from correction_factor.correction_factor_2 import correction_factor_given_radex_table

# Default values come from Goldsmith/Bergin/Lis '97, esp. the c18o abundance
def make_goldsmith_figure_9(save=True, column=1e16, abundance=1.7e-7, print_timing=False):
    """ Generates my version of Fig. 9 from Goldsmith+'97 """

    start = time.time()

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

    end = time.time()
    if print_timing:
        print end-start

    return fig


def make_goldsmith_figure_10(save=True, print_timing=False):
    """ Generates my version of Fig. 10 from Goldsmith+'97 """

    start = time.time()

    # Here we have to specify 2 of 3: abundance, temperature, density.
    R = pyradex.Radex(species='c18o', abundance=1.7e-7, column=1e14, temperature=50)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_array = np.logspace(1, 7, 100)
    temperature_list = [10, 30, 50]

    f_c_array_10K = np.zeros_like(density_array)
    f_c_array_30K = np.zeros_like(density_array)
    f_c_array_50K = np.zeros_like(density_array)

    f_c_array_list = [f_c_array_10K, f_c_array_30K, f_c_array_50K]

    for j, temp in enumerate(temperature_list):

        R.temperature = temp

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

    end = time.time()
    if print_timing:
        print end-start

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


def make_plume_figure_3(save=True, print_timing=False, abundance=5e-7, n_points=20):
    """ Generates my version of Fig. 3 from Plume+'12 """

    # Plume '12 claims to use an abundance of 2e-7 (i.e. one C18O for every 500 x 10^4 H2 atoms)
    # but I find that 5e-7 matches the appearance much more closely.

    start = time.time()

    # Figure 3's caption says that Ju=10 is included, but the text omits it. 
    # I have had better luck reproducing the plot when I omit Ju=10,
    # so I infer that the caption is incorrect.
    J_upper_list = [2,3,5,7,8,9,11,15]
    J_lower_list = [x-1 for x in J_upper_list]

    # Here we have to specify 2 of 3: abundance, temperature, density.
    R = pyradex.Radex(species='c18o', abundance=abundance, column=5e16, temperature=50, escapeProbGeom='lvg')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_list = [1e5, 1e6, 1e7]
    temperature_array = np.linspace(10, 400, n_points)

    f_c_array_n5 = np.zeros_like(temperature_array)
    f_c_array_n6 = np.zeros_like(temperature_array)
    f_c_array_n7 = np.zeros_like(temperature_array)

    f_c_array_list = [f_c_array_n5, f_c_array_n6, f_c_array_n7]

    for j, density in enumerate(density_list):

        R.density = density

        for i, temp in enumerate(temperature_array):

            new_T = R(temperature=temp)

            # use the correction_factor code here!
            f_c = correction_factor_given_radex_table(new_T, J_lower_list, included_levels=18)
            # They only use the lowest 18 rotational levels, 
            # even though the plot's appearance changes when 40 are included.

            f_c_array_list[j][i] = f_c
            
    ax.plot(temperature_array, f_c_array_n7, 'b-', lw=1.5, label=r'$n = 10^7\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n6, 'g--', lw=1.5, label=r'$n = 10^6\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n5, 'r:', lw=2, label=r'$n = 10^5\ \rm{cm}^{-3}$')

    ax.legend(frameon=False, loc='upper center', fontsize=18)

    ax.set_xlim(0, 400)
    ax.set_ylim(1.6, 2)
    ax.set_xticks([0,100,200,300,400])
    ax.set_yticks([1.6,1.7,1.8,1.9,2])

    fontdict = {'family':'serif', 'fontsize':16}

    ax.set_ylabel("Correction Factor", fontdict=fontdict)
    ax.set_xlabel("Temp (K)", fontdict=fontdict)

    ax.minorticks_on()

    if save:
        fig.savefig("plume12_fig_3_reproduction.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print end-start

    return fig    