"""
A function to build an astropy table out of the myradex output.

"""

from __future__ import division

from astropy.table import Table

# The following line is printed when myradex initializes, so I'm gonna trust it.
# Column names of the output: iup ilow Eup freq lam Tex tau Tr fup flow flux_K flux beta Jnu gup glow Aul Bul Blu Jback flux_dens

string_of_colnames = "iup ilow Eup freq lam Tex tau Tr fup flow flux_K flux beta Jnu gup glow Aul Bul Blu Jback flux_dens"
list_of_colnames = string_of_colnames.split(" ")

def table_from_myradex_datatransitions(data_transitions):

    if data_transitions.shape[0] != len(list_of_colnames):
        raise ValueError("`data_transitions` has incorrect number of columns. ")

    output_table = Table(data=data_transitions.T, names=list_of_colnames)

    return output_table