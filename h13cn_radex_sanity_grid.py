"""
This module is developed in response to the fact that some things aren't lining up nicely.

Does RADEX treat h13cn inconsistently?

""" 

from __future__ import division
import time

import numpy as np
import matplotlib.pyplot as plt
import pyradex

# This import will have to change if I refactor correction_factor much.
from correction_factor.correction_factor_2 import correction_factor_given_radex_table

h13cn_J_upper_list = [1, 3, 4, 6, 7, 8, 9, 10]
h13cn_J_lower_list = [x-1 for x in h13cn_J_upper_list]

# ok so let's g

n_temperatures = 8
min_temp = 20
max_temp = 400

m_logdensities = 8
min_logdens = 4
max_logdens = 9

temperatures = np.linspace(min_temp, max_temp, n_temperatures)
logdensities = np.linspace(min_logdens, max_logdens, m_logdensities)

sanity_grid = np.zeros((n_temperatures, m_logdensities))
iters_grid = np.zeros((n_temperatures, m_logdensities))

# Initialize radex object
initial_column=1e16 
initial_abundance=7e-8
R = pyradex.Radex(species='h13cn@xpol', abundance=initial_abundance, temperature=20, 
                  column=initial_column)
R.maxiter = 1000

abundance_list = [2e-9, 2e-8, 7e-8]

for abundance in abundance_list:

    R.abundance = abundance

    for i, temp in enumerate(temperatures):

        R.temperature = temp

        for j, logdens in enumerate(logdensities):

            R.density = 10**logdens

            # calculate the T thing

            if j == 0:
                n_iters = R.run_radex(reuse_last=False, reload_molfile=True, silent=False)
            else:
                n_iters = R.run_radex()

            T = R.get_table()

            print "Tkin = {0}, logdens = {1}, column = {2}".format(R.temperature, 
                                                                   np.log10(R.density['H2'].value), 
                                                                   R.column)
            print "Tex of 1-0 line: {0}".format(T['Tex'][0])
            print "n_iters={1}. Sum of lower level pops: {0}".format(np.sum(T['lowerlevelpop']), 
                                                                     n_iters)
            print ""

            negative_tex_count = len(T[T['Tex'] < 0])

            sanity_grid[i,j] = negative_tex_count
            iters_grid[i,j] = n_iters


            # if (min(T['Tex']) < 0) and (np.sum(T['lowerlevelpop']) < 0.99):
            #     # this will FAIL on last iteration if I mixed up my axis order
            #     sanity_grid[i,j] = -2
            # elif min(T['Tex']) < 0:
            #     sanity_grid[i,j] = -1
            # elif np.sum(T['lowerlevelpop']) < 0.99:
            #     sanity_grid[i,j] = 0
            # else:
            #     sanity_grid[i,j] = 1


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(sanity_grid, origin='lower', extent=(min_logdens, max_logdens, min_temp, max_temp),
              interpolation='nearest', aspect=(max_logdens-min_logdens)/(max_temp-min_temp), 
              cmap='Reds')

    fontdict = {'family':'serif', 'fontsize':16}

    ax.set_xlabel("log (H$_2$ density/cm$^{-3}$)", fontdict=fontdict)
    ax.set_ylabel("Temp (K)", fontdict=fontdict)
    ax.set_title("H$^{13}$CN Column density: %.1e, abundance=%.1e" % (R.column.value, R.abundance))

    plt.savefig("X{0:1e}.png".format(abundance))

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.imshow(iters_grid, origin='lower', extent=(min_logdens, max_logdens, min_temp, max_temp),
              interpolation='nearest', aspect=(max_logdens-min_logdens)/(max_temp-min_temp), 
              cmap='Blues')

    fontdict = {'family':'serif', 'fontsize':16}

    ax2.set_xlabel("log (H$_2$ density/cm$^{-3}$)", fontdict=fontdict)
    ax2.set_ylabel("Temp (K)", fontdict=fontdict)
    ax2.set_title("H$^{13}$CN Column density: %.1e, abundance=%.1e" % (R.column.value, R.abundance))

    plt.show()

