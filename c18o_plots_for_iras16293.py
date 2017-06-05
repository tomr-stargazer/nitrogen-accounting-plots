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
# from c18o_myradex_interface import myradex_c18o

def c18o_cf_iras16293(save=True, print_timing=False, abundance=5e-7, n_points=20):
    """ Generates my version of Fig. 3 from Plume+'12 """

    # Plume '12 claims to use an abundance of 2e-7 (i.e. one C18O for every 500 x 10^4 H2 atoms)
    # but I find that 5e-7 matches the appearance much more closely.

    start = time.time()

    # Figure 3's caption says that Ju=10 is included, but the text omits it. 
    # I have had better luck reproducing the plot when I omit Ju=10,
    # so I infer that the caption is incorrect.
    J_upper_list = [2,3,5,6,7,8,9,10]
    J_lower_list = [x-1 for x in J_upper_list]

    # Here we have to specify 2 of 3: abundance, temperature, density.
    R = pyradex.Radex(species='c18o', abundance=abundance, column=5e16, temperature=50, escapeProbGeom='lvg')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_list = [1e5, 1e6, 1e7, 1e8]
    temperature_array = np.linspace(10, 400, n_points)

    f_c_array_n5 = np.zeros_like(temperature_array)
    f_c_array_n6 = np.zeros_like(temperature_array)
    f_c_array_n7 = np.zeros_like(temperature_array)
    f_c_array_n8 = np.zeros_like(temperature_array)

    f_c_array_list = [f_c_array_n5, f_c_array_n6, f_c_array_n7, f_c_array_n8]

    for j, density in enumerate(density_list):

        R.density = density

        for i, temp in enumerate(temperature_array):

            new_T = R(temperature=temp)

            # use the correction_factor code here!
            f_c = correction_factor_given_radex_table(new_T, J_lower_list)

            f_c_array_list[j][i] = f_c

    ax.plot(temperature_array, f_c_array_n8, '-.', color='purple', lw=2, label=r'$n = 10^8\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n7, 'b-', lw=1.5, label=r'$n = 10^7\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n6, 'g--', lw=1.5, label=r'$n = 10^6\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n5, 'r:', lw=2, label=r'$n = 10^5\ \rm{cm}^{-3}$')

    ax.legend(frameon=False, loc='upper center', fontsize=18)

    ax.set_xlim(0, 400)
    ax.set_ylim(1.2, 2.4)
    ax.set_xticks([0,100,200,300,400])
    # ax.set_yticks([1.6,1.7,1.8,1.9,2])

    fontdict = {'family':'serif', 'fontsize':16}

    ax.set_ylabel("Correction Factor", fontdict=fontdict)
    ax.set_xlabel("Temp (K)", fontdict=fontdict)

    ax.minorticks_on()

    if save:
        fig.savefig("c18o_iras16293.pdf")
        np.savetxt("c18o_fc_v_T_myr_temperature_array.txt", temperature_array)
        np.savetxt("c18o_fc_v_T_myr_fc_array_n5.txt", f_c_array_n5)
        np.savetxt("c18o_fc_v_T_myr_fc_array_n6.txt", f_c_array_n6)
        np.savetxt("c18o_fc_v_T_myr_fc_array_n7.txt", f_c_array_n7)
        np.savetxt("c18o_fc_v_T_myr_fc_array_n8.txt", f_c_array_n8)

    plt.show()

    end = time.time()
    if print_timing:
        print end-start

    return fig    