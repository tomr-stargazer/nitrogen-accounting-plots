"""
Here is some functionality around making sure that the RADEX output table results are sane.

"""

from __future__ import division

import numpy as np


def check_table_for_sanity():
    pass

def Tex_inspection():
    """ Ensure all Tex values are positive. """

    pass

def population_sum_inspection(radex_table, tolerance=0.01):
    """ Ensure that the fractional population of levels sums to ~1. """

    if np.sum(radex_table['lowerlevelpop']) < 1-tolerance:
        return False
    else:
        return True