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
from collections import defaultdict
import numpy as np
import matplotlib as mpl
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


def combine_XCLASS_MADEX_preferring_MADEX():

    madex = load_table_8_MADEX()
    xclass = load_table_8_XCLASS()

    # This loop modifies table `xclass` in-place.
    for i in range(len(madex)):

        xclass[xclass['Molecule'] == madex['Molecule'][i]] = madex[i]

    return xclass


def load_table_7_Zernickel():

    table7_zernickel = astropy.table.Table.read("zernickel_table7_just_the_table.html", header_start=0, data_start=4)

    return table7_zernickel


def transform_table_7_Zernickel(table):

    # so basically, we wanna ditch the other columns and get out abundances (relative to N_H2 like a normal person).

    Xz_values = np.array([float(x.replace("(","e").replace("–","-").rstrip(")")) for x in table['NGC 6334I']])
    Xs = astropy.table.Column(2*Xz_values, name="Abundance (X/H2)")

    new_table = astropy.table.Table([table['Molecule'], Xs])

    return new_table


def latex_molecule_name(molecule_name):
    """ Makes the string latex-friendly. """
    latex_name = r"$\rm{"+molecule_name+"}$"
    return latex_name


def latex_molecule_name_zern(molecule_name):
    """ Makes the Zernickel name string latex-friendly. """
    integers = [str(x) for x in range(9)]
    for integer in integers:
        molecule_name = molecule_name.replace(integer, '_'+integer)

    latex_name = r"$\rm{"+molecule_name+"}$"
    return latex_name


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
orion_KL_ylims = (0.0003, 1.5)
orion_KL_yticks = ((3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 0.3, 0.6, 1), 
                   (r"0.03\%", r"0.1\%", r"0.3\%", r"1\%", r"3\%", r"10\%", r"30\%", r"60\%",r"100\%"))
ngc6334_colors = ['tab:red', 'tab:purple', 'tab:olive', 'tab:green', 'tab:pink', 'tab:gray']
ngc6334_savename = "zern_new_plot.png"
ngc6334_ylims = (0.0003, 1.5)


def plot_organic_nitrogen_fraction_w_colors_and_errors(
    table, abundance_colname='Hot Core', ylims=orion_KL_ylims, yticks=orion_KL_yticks, colors=orion_KL_colors, 
    savename=orion_KL_savename, dpi=400, latex_molname_fn=latex_molecule_name):
    """ starting basic then adding bells-and-whistles like sorting. """

    fig = plt.figure(figsize=(2.5, 2.5), dpi=dpi)
    ax = fig.add_subplot(111)

    organic_N_table = table[['N' in x['Molecule'] and 'C' in x['Molecule'] and ~np.isnan(x[abundance_colname]) for x in table]]
    organic_N_table.sort(abundance_colname)
    organic_N_table.reverse()

    naive_total_organic_N = np.sum(organic_N_table[abundance_colname])
    print("Total organic N: {0:.3e}".format(naive_total_organic_N))

    lx_molnames = [latex_molname_fn(x) for x in organic_N_table['Molecule']]

    # colors = ['tab:red', 'tab:purple', 'tab:green', 'tab:olive', 'tab:cyan', 'tab:pink', 'tab:brown', 'tab:orange']
    for i in range(len(organic_N_table)):

        color = colors[i]

        total_organic_N_others = naive_total_organic_N - organic_N_table[abundance_colname][i]
        fraction = organic_N_table[abundance_colname][i]/naive_total_organic_N
        abundance_high = organic_N_table[abundance_colname][i]*1.35
        fraction_high = abundance_high/(total_organic_N_others+abundance_high)
        abundance_low = organic_N_table[abundance_colname][i]*0.65
        fraction_low = abundance_low/(total_organic_N_others+abundance_low)

        if True: #organic_N_table['Molecule'][i] == 'HCN':
            print("{0}:".format(organic_N_table['Molecule'][i]))
            print(fraction, fraction_high, fraction_low)
            print(organic_N_table[abundance_colname][i], abundance_high, abundance_low)

        ax.errorbar(i, fraction, yerr=np.array([fraction-fraction_low, fraction_high-fraction]).reshape(2,1), 
                    fmt='o', ms=3.5, color=color, ecolor=color, capsize=3, elinewidth=1.5, capthick=1.5, zorder=10)

        t = ax.text(i+0.25, fraction*1.75, lx_molnames[i], bbox={'facecolor':'white', 'edgecolor':'none', 'pad':0.5},
                    fontdict={'family':'serif', 'weight':'bold', 'color':color, 'rotation':30}, zorder=5)


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


def generate_molecule_to_color_dict(table, input_colors, abundance_colname='Hot Core'):

    organic_N_table = table[['N' in x['Molecule'] and 'C' in x['Molecule'] and ~np.isnan(x[abundance_colname]) for x in table]]
    organic_N_table.sort(abundance_colname)
    organic_N_table.reverse()

    mol_color_dict = defaultdict(lambda : 'tab:gray')

    for i in range(len(organic_N_table)):

        key = organic_N_table['Molecule'][i].replace('_', '')

        mol_color_dict[key] = input_colors[i]

    return mol_color_dict


def colorlist_from_molcolordict(table, mol_color_dict, abundance_colname='Abundance (X/H2)'):

    organic_N_table = table[['N' in x['Molecule'] and 'C' in x['Molecule'] and ~np.isnan(x[abundance_colname]) for x in table]]
    organic_N_table.sort(abundance_colname)
    organic_N_table.reverse()

    color_list = [mol_color_dict[x.replace('_', '')] for x in organic_N_table['Molecule']]

    return color_list


def plot_organic_nitrogen_fraction_w_colors_and_errors_zernickel(savename=ngc6334_savename, colors=ngc6334_colors, **kwargs):

    table = transform_table_7_Zernickel(load_table_7_Zernickel())

    return plot_organic_nitrogen_fraction_w_colors_and_errors(
         table, abundance_colname='Abundance (X/H2)', colors=colors, ylims=ngc6334_ylims, savename=savename, **kwargs)


def get_N_molecules(table):
    """ Assumes table has a 'Molecule' column. """

    table8[['N' in x for x in table8['Molecule']]]

    molecules_unsorted = np.array([molecule for molecule, value in items])
    values_unsorted = np.array([value for molecule, value in items])
    sorted_indices = np.argsort(values_unsorted)[::-1]

    molecules_sorted_by_max_value = molecules_unsorted[sorted_indices]





def pie_organic_nitrogen_fraction(table, abundance_colname='Hot Core', colors=None, latex_molname_fn=latex_molecule_name,
                                  dx_array=None, dy_array=None, rotation_array=None):
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
    labels = [latex_molname_fn(x) for x in organic_N_table['Molecule']]

    explode_array = np.ones_like(fractions) * 0.1
    explode_array[0] = 0

    slices, texts = ax.pie(fractions, labels=labels, colors=colors, startangle=45, explode=explode_array,
           wedgeprops={'linewidth':0.5, 'ec':'k'}, counterclock=False, textprops={'family':'serif', 'weight':'bold'})

    # I will *probably* have to implement my own labels rather than fighting to get these ones to work like this.
    # OR... I'll have to modify the text objects post-facto!

    for i, t in enumerate(texts):

        t.set_color(colors[i])

        if dx_array is not None and dy_array is not None:
            pass
        if rotation_array is not None:
            pass

    return fig


if __name__ == "__main__":

    import os

    figure_dir = os.path.expanduser("~/Documents/Academia/Articles/Nitrogen_Paper/in_progress_graphics/")

    xm = combine_XCLASS_MADEX_preferring_MADEX()
    colors = ['C{0}'.format(x) for x in range(10)]
    orion_KL_colors = colors # ['tab:red', 'tab:blue', 'tab:green', 'tab:purple', 'tab:cyan', 'tab:olive', 'tab:pink', 'tab:brown', 'tab:orange']
    # too bad colors isn't a dict

    plot_organic_nitrogen_fraction_w_colors_and_errors(xm, 
                                                       savename=os.path.join(figure_dir, "Orion_KL_accounting_logfraction.pdf"), 
                                                       colors=orion_KL_colors, dpi=100)

    okl_pie = pie_organic_nitrogen_fraction(xm, colors=orion_KL_colors)
    okl_pie.savefig(os.path.join(figure_dir, "Orion_KL_accounting_pie.pdf"), bbox_inches='tight')


    zern_table = transform_table_7_Zernickel(load_table_7_Zernickel())
    mol_color_dict = generate_molecule_to_color_dict(xm, orion_KL_colors)
    zern_colors = colorlist_from_molcolordict(zern_table, mol_color_dict)

    plot_organic_nitrogen_fraction_w_colors_and_errors_zernickel(savename=os.path.join(figure_dir, "NGC6334I_accounting_logfraction.pdf"), 
                                                                 colors=zern_colors, dpi=100, latex_molname_fn=latex_molecule_name_zern)

    zern_pie = pie_organic_nitrogen_fraction(zern_table, colors=zern_colors, abundance_colname='Abundance (X/H2)', latex_molname_fn=latex_molecule_name_zern)
    zern_pie.savefig(os.path.join(figure_dir, "NGC6334I_accounting_pie.pdf"), bbox_inches='tight')

    plt.show()


    # do 2 plots for each source
