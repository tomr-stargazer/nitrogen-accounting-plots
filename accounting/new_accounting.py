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

    organic_N_table = table[['N' in x and 'C' in x for x in table['Molecule']]]
    organic_N_table.sort("Hot Core")
    organic_N_table.reverse()

    inorganic_N_table = table[['N' in x and 'C' not in x for x in table['Molecule']]]
    inorganic_N_table.sort("Hot Core")
    inorganic_N_table.reverse()

    N_table = astropy.table.vstack([organic_N_table, inorganic_N_table])

    ax.plot(N_table['Hot Core'], 'ro')
    ax.set_xticks(np.arange(len(N_table)))
    ax.set_xticklabels(N_table['Molecule'], rotation=90)

    ax.semilogy()

    xs = len(organic_N_table)-0.5 * np.ones(2)
    ys = [1e-10, 1e-2]
    ax.plot(xs, ys, 'k--', scaley=False, scalex=False)

    ax.set_ylabel("Abundance (X/H$_2$)")

    return fig


def plot_organic_nitrogen_fraction(table):
    """ starting basic then adding bells-and-whistles like sorting. """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    organic_N_table = table[['N' in x and 'C' in x for x in table['Molecule']]]
    organic_N_table.sort("Hot Core")
    organic_N_table.reverse()

    total_organic_N = np.nansum(organic_N_table['Hot Core'])

    ax.plot(organic_N_table['Hot Core']/total_organic_N, 'ro')
    ax.set_xticks(np.arange(len(organic_N_table)))
    ax.set_xticklabels(organic_N_table['Molecule'], rotation=90)

    ax.set_ylabel("Organic N fraction")

    ax.semilogy()
    ax.set_ylim(0.005, 1.5)
    plt.yticks(
        (1e-2, 3e-2, 1e-1, 0.3, 0.6, 1), 
        (r"1%", r"3%", r"10%", r"30%", r"60%",r"100%"), fontsize=8)

    return fig


def get_N_molecules(table):
    """ Assumes table has a 'Molecule' column. """

    table8[['N' in x for x in table8['Molecule']]]

    molecules_unsorted = np.array([molecule for molecule, value in items])
    values_unsorted = np.array([value for molecule, value in items])
    sorted_indices = np.argsort(values_unsorted)[::-1]

    molecules_sorted_by_max_value = molecules_unsorted[sorted_indices]


