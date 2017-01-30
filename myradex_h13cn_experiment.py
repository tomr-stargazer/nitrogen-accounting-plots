"""
Experimenting with myradex.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import wrapper_my_radex
wrapper = wrapper_my_radex.myradex_wrapper

n_levels, n_item, n_transitions = wrapper.config_basic(
    '/Users/tsrice/Documents/Academia/Misc_Software/Radex/data/',
    'h13cn@xpol.dat', 2.7, True)

params = {'tkin': 1e3,
          'dv_cgs': 1e5,
          'dens_x_cgs': 1e6,
          'ncol_x_cgs': 1e15,
          'h2_density_cgs': 1e9,
          'hi_density_cgs': 1e1,
          'oh2_density_cgs': 0.0,
          'ph2_density_cgs': 0.0,
          'hii_density_cgs': 0.0,
          'electron_density_cgs': 1e6,
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
