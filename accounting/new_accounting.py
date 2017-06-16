"""
Ok, so the results I showed in Chile using the following scripts:
- N_barcharts.py
- astronomical_NSi_ratios_poster.py
were totally actually wrong. 

Specifically, I was trying to back abundances out of Table 4 of Crockett+14b, 
but totally ignored (a) the discussion in Sections 4 and 5 of that paper 
describing how to measure abundances generally and molecule-by-molecule and 
(b) HEY TOM TABLE 8 HAS THE HECKIN ABUNDANCES RIGHT THERE???

So anyway, this code is for doing it right.

"""

import numpy as np
import matplotlib.pyplot as plt

import astropy.table

def load_table_8_XCLASS():

    table8_xclass = astropy.table.Table.read(
        "apj493753t8_ascii_copy.txt", format='ascii.basic', delimiter='\t', 
        guess=False, header_start=2, data_start=4, data_end=37)

    return table8_xclass


def load_table_8_MADEX():

    table8_madex = astropy.table.Table.read(
        "apj493753t8_ascii_copy.txt", format='ascii.basic', delimiter='\t', 
        guess=False, header_start=2, data_start=39, data_end=44)

    return table8_madex


def plot_molecular_abundances(table):
    """ starting basic then adding bells-and-whistles like sorting. """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    N_table = table[['N' in x for x in table['Molecule']]]
    N_table.sort("Hot Core")
    N_table.reverse()

    ax.plot(N_table['Hot Core'], 'ro')
    ax.set_xticks(np.arange(len(N_table)))
    ax.set_xticklabels(N_table['Molecule'], rotation=90)

    ax.semilogy()

    return fig


def get_N_molecules(table):
    """ Assumes table has a 'Molecule' column. """

    table8[['N' in x for x in table8['Molecule']]]

    molecules_unsorted = np.array([molecule for molecule, value in items])
    values_unsorted = np.array([value for molecule, value in items])
    sorted_indices = np.argsort(values_unsorted)[::-1]

    molecules_sorted_by_max_value = molecules_unsorted[sorted_indices]


