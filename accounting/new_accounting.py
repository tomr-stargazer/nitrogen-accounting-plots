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

import pdb
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


def load_table_7_Zernickel():

    table7_zernickel = astropy.table.Table.read("zernickel_table7_just_the_table.html", header_start=0, data_start=4)

    return table7_zernickel


def transform_table_7_Zernickel(table):

    # so basically, we wanna ditch the other columns and get out abundances (relative to N_H2 like a normal person).

    Xz_values = np.array([float(x.replace("(","e").replace("–","-").rstrip(")")) for x in table['NGC 6334I']])
    Xs = astropy.table.Column(2*Xz_values, name="Abundance (X/H2)")

    new_table = astropy.table.Table([table['Molecule'], Xs])

    return new_table


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
    ax.set_xticks([])

    ax.semilogy()
    ax.set_ylim(0.005, 1.5)
    plt.yticks(
        (1e-2, 3e-2, 1e-1, 0.3, 0.6, 1), 
        (r"1%", r"3%", r"10%", r"30%", r"60%",r"100%"), fontsize=8)

    return fig


orion_KL_colors = ['tab:red', 'tab:purple', 'tab:green', 'tab:olive', 'tab:cyan', 'tab:pink', 'tab:brown', 'tab:orange']
orion_KL_savename = "new_nitrogen_fraction_plot_for_poster.png"
orion_KL_ylims = (0.003, 1.5)
orion_KL_yticks = ((0.003, 1e-2, 3e-2, 1e-1, 0.3, 0.6, 1), 
                   (r"0.3%", r"1%", r"3%", r"10%", r"30%", r"60%",r"100%"))
ngc6334_colors = ['tab:red', 'tab:purple', 'tab:olive', 'tab:green', 'tab:pink', 'tab:gray']
ngc6334_savename = "zern_new_plot.png"
ngc6334_ylims = (0.003, 1.5)


def plot_organic_nitrogen_fraction_w_colors_and_errors(
    table, abundance_colname='Hot Core', ylims=orion_KL_ylims, yticks=orion_KL_yticks, colors=orion_KL_colors, 
    savename=orion_KL_savename):
    """ starting basic then adding bells-and-whistles like sorting. """

    fig = plt.figure(figsize=(2.5, 2.5), dpi=400)
    ax = fig.add_subplot(111)

    organic_N_table = table[['N' in x['Molecule'] and 'C' in x['Molecule'] and ~np.isnan(x[abundance_colname]) for x in table]]
    organic_N_table.sort(abundance_colname)
    organic_N_table.reverse()

    naive_total_organic_N = np.sum(organic_N_table[abundance_colname])

    # colors = ['tab:red', 'tab:purple', 'tab:green', 'tab:olive', 'tab:cyan', 'tab:pink', 'tab:brown', 'tab:orange']
    for i, color in enumerate(colors):

        total_organic_N_others = naive_total_organic_N - organic_N_table[abundance_colname][i]
        fraction = organic_N_table[abundance_colname][i]/naive_total_organic_N
        abundance_high = organic_N_table[abundance_colname][i]*1.4
        fraction_high = abundance_high/(total_organic_N_others+abundance_high)
        abundance_low = organic_N_table[abundance_colname][i]*0.6
        fraction_low = abundance_low/(total_organic_N_others+abundance_low)

        if organic_N_table['Molecule'][i] == 'HCN':
            print(fraction, fraction_high, fraction_low)
            print(organic_N_table[abundance_colname][i], abundance_high, abundance_low)

        ax.errorbar(i, fraction, yerr=np.array([fraction-fraction_low, fraction_high-fraction]).reshape(2,1), 
                    fmt='o', ms=3.5, color=color, ecolor=color, capsize=3, elinewidth=1.5, capthick=1.5)


    # ax.plot(organic_N_table[abundance_colname]/total_organic_N, 'ro')
    ax.set_xticks(np.arange(len(organic_N_table)))
    ax.set_xticklabels(organic_N_table['Molecule'], rotation=90)

    ax.set_ylabel("Organic N fraction")
    ax.set_xticks([])

    ax.semilogy()
    ax.set_ylim(ylims)
    plt.gca().yaxis.grid()
    plt.yticks(*yticks, fontsize=8)

    if savename is not None:
        plt.savefig(savename, bbox_inches='tight')

    return fig


def plot_organic_nitrogen_fraction_w_colors_and_errors_zernickel(savename=ngc6334_savename):

    table = transform_table_7_Zernickel(load_table_7_Zernickel())

    return plot_organic_nitrogen_fraction_w_colors_and_errors(
         table, abundance_colname='Abundance (X/H2)', colors=ngc6334_colors, ylims=ngc6334_ylims, savename=savename)


def get_N_molecules(table):
    """ Assumes table has a 'Molecule' column. """

    table8[['N' in x for x in table8['Molecule']]]

    molecules_unsorted = np.array([molecule for molecule, value in items])
    values_unsorted = np.array([value for molecule, value in items])
    sorted_indices = np.argsort(values_unsorted)[::-1]

    molecules_sorted_by_max_value = molecules_unsorted[sorted_indices]


def latex_molecule_name(molecule_name):
    """ Makes the string latex-friendly. """
    latex_name = r"$\rm{"+molecule_name+"}$"
    return latex_name


def pie_organic_nitrogen_fraction(table, abundance_colname='Hot Core', colors=None):
    """ starting basic then adding bells-and-whistles like sorting. """

    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)

    organic_N_table = table[['N' in x['Molecule'] and 'C' in x['Molecule'] and ~np.isnan(x[abundance_colname]) for x in table]]
    organic_N_table.sort(abundance_colname)
    organic_N_table.reverse()    

    total_organic_N = np.nansum(organic_N_table[abundance_colname])

    # labels = "HCN jk", "CH3CN nope", "something", "yes"
    # fractions = [0.45, 0.2, 0.2, 0.15]

    fractions = organic_N_table[abundance_colname]/total_organic_N
    labels = [latex_molecule_name(x) for x in organic_N_table['Molecule']]

    ax.pie(fractions, labels=labels, colors=colors)

    # ax.plot(organic_N_table['Hot Core']/total_organic_N, 'ro')
    # ax.set_xticks(np.arange(len(organic_N_table)))
    # ax.set_xticklabels(organic_N_table['Molecule'], rotation=90)

    # ax.set_ylabel("Organic N fraction")
    # ax.set_xticks([])

    # ax.semilogy()
    # ax.set_ylim(0.005, 1.5)
    # plt.yticks(
    #     (1e-2, 3e-2, 1e-1, 0.3, 0.6, 1), 
    #     (r"1%", r"3%", r"10%", r"30%", r"60%",r"100%"), fontsize=8)

    return fig
