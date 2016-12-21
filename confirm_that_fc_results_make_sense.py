"""
This module is developed in response to the fact that some things aren't lining up nicely.

In some of my plots, it seems to make a big difference whether we go "density-first" or 
"temperature-first" when iterating through RADEX computations that are used to get f_c numbers.

The specific observation that seems to be upsetting:

when the f_c of H13CN is computed (for J_u list = [1, 3, 4, 6, 7, 8, 9, 10]), at n=10^6 
and T=50, 100, 200 K, I get (starkly) different numbers depending on whether I am iterating
over temperature first or density first. I trust the "temperature first" results a little better.

I don't understand RADEX very well.

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

def crude_attempt_at_matching():
    """ 
    Here we are trying to... convince correction_factor and RADEX to be 
    consistent across different kinds of runs!

    Step 1: do temperature then density
    Step 2: do density then temperature

    Step 3: make sure they are consistent dude
    Step 4: do them in reverse order and try again

    I have a vague sense that (because stuff is directly interfacing with 
    underlying FORTRAN code) 

    """

    # Here we have to specify 2 of 3: abundance, temperature, density.
    R_step1 = pyradex.Radex(species='h13cn@xpol', abundance=1e-8, temperature=100, 
                            column=1e16, escapeProbGeom='lvg')

    density = 1e6
    temperature_list = [50, 100, 200]

    f_c_dict = {}

    for j, temp in enumerate(temperature_list):

        R_step1.temperature = temp

        new_T = R_step1(collider_densities={'H2':density})

        f_c = correction_factor_given_radex_table(new_T, h13cn_J_lower_list)

        f_c_dict[temp] = f_c

    print f_c_dict

    del R_step1
                
    # Step 2
    R_step2 = pyradex.Radex(species='h13cn@xpol', abundance=1e-8, column=1e16, temperature=50, escapeProbGeom='lvg')

    f_c_dict_2 = {}

    R_step2.density = density

    for i, temp in enumerate(temperature_list):

        new_T = R_step2(temperature=temp)

        # use the correction_factor code here!
        f_c = correction_factor_given_radex_table(new_T, h13cn_J_lower_list)

        f_c_dict_2[temp] = f_c

    print f_c_dict_2


    return 