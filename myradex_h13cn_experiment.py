"""
Experimenting with myradex.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import wrapper_my_radex
wrapper = wrapper_my_radex.myradex_wrapper

from correction_factor.correction_factor_2 import correction_factor_given_myradex_output

import h13cn_myradex_interface

n_levels, n_item, n_transitions = wrapper.config_basic(
    '/Users/tsrice/Documents/Academia/Misc_Software/Radex/data/',
    'h13cn@xpol.dat', 2.7, True)

params = {'tkin': 400,
          'dv_cgs': 1e5, # What is this? oh, it's delta v in cm/scp
          'dens_x_cgs': 1e2, # I think this is not a relevant parameter?
          'ncol_x_cgs': 5e16,
          'h2_density_cgs': 1e6,
          'hi_density_cgs': 0.0,
          'oh2_density_cgs': 0.0,
          'ph2_density_cgs': 0.0,
          'hii_density_cgs': 0.0,
          'electron_density_cgs': 0.0,
          'n_levels': n_levels,
          'n_item': n_item,
          'n_transitions': n_transitions,
          'geotype':'lvg'}

energies,f_occupations,data_transitions,cooling_rate = \
    wrapper.run_one_params(**params)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(energies, f_occupations, marker='o')
ax.set_xlabel('E (K)')
ax.set_ylabel('Occupation')
ax.set_ylim((1e-15,1e0))

observed_lines = [1, 3, 4]

print correction_factor_given_myradex_output(f_occupations, observed_lines)


# now let's do something fancier

temp_array = [10, 80, 200, 400]

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

for T in temp_array:

    loop_params = params.copy()
    loop_params['tkin'] = T

    energies,f_occupations,data_transitions,cooling_rate = \
        wrapper.run_one_params(**loop_params)

    ax2.plot(energies, f_occupations, marker='o')

    print "fcorr at T={0}:".format(T)
    print correction_factor_given_myradex_output(f_occupations, observed_lines)


ax2.set_xlabel('E (K)')
ax2.set_ylabel('Occupation')


plt.show()
