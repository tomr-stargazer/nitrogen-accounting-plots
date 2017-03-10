"""
We're gonna try and make a bar-chart oriented dealio here.

"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import astropy.table


def make_preliminary_figure(table):
    """
    A prototype. Playing with summing and sorting. Maybe with display, too.

    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    molecule_list = list(set(table['Molecule']))
    unsorted_abundances = np.array([np.sum(table['N_tot'][table['Molecule']==molecule]) for molecule in molecule_list])
    sorted_indices = np.argsort(unsorted_abundances)[::-1]
    sorted_abundances = unsorted_abundances[sorted_indices]
    sorted_molecule_list = np.array(molecule_list)[sorted_indices]

    plt.plot(sorted_abundances, 'ko')

    # xticklabels = []
    # for i, molecule in enumerate(molecule_set):
    #     ax.plot([i], [np.sum(table['N_tot'][table['Molecule']==molecule])], 'ko')
    #     xticklabels.append(molecule)

    ax.set_xticks(np.arange(len(molecule_list)))
    ax.set_xticklabels(sorted_molecule_list, rotation=90)
    plt.show()
    return fig


def main_species(molecule_name):
    """ Uses simple rules to isolate main species """
    name_sans_vibrational_state = molecule_name.split(',')[0]
    name_sans_13C = name_sans_vibrational_state.replace("^13C", 'C')
    name_sans_15N = name_sans_13C.replace("^15N", 'N')

    return name_sans_15N


def make_coldens_range_dict(table):
    """ result: (main isotopologue -> tuple of coldens values) dict """

    ratio_12C_13C = 69
    ratio_14N_15N = 388

    molecule_dict = {}

    for row in table:
        # is this an isotopologue or a main species?
        if '^13C' in row['Molecule']:
            N_tot_main = row['N_tot'] * ratio_12C_13C
        elif '^15N' in row['Molecule']:
            N_tot_main = row['N_tot'] * ratio_14N_15N
        else:
            N_tot_main = row['N_tot']

        try:
            molecule_dict[main_species(row['Molecule'])].append(N_tot_main)
        except KeyError:
            molecule_dict[main_species(row['Molecule'])] = [N_tot_main]
        
    return molecule_dict



def make_errorbar_figure(table):

    molecules_in_order = ['HCN', 'NH3',]

    h2_column_density = 3.1e23

    fig = plt.figure()
    ax = add_subplot(111)

    ax.set_ylabel("Abundance")

    return fig


