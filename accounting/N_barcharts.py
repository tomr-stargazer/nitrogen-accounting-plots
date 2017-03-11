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


def latex_molecule_name(molecule_name):
    """ Makes the string latex-friendly. """
    latex_name = r"$\rm{"+molecule_name+"}$"
    return latex_name


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


def make_errorbar_figure(molecule_dict, ylabel="Column density (cm$^{-2}$)"):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    max_value_per_molecule_dict = {x: np.max(y) for x, y in molecule_dict.items()}
    items = max_value_per_molecule_dict.items()
    molecules_unsorted = np.array([molecule for molecule, value in items])
    values_unsorted = np.array([value for molecule, value in items])
    sorted_indices = np.argsort(values_unsorted)[::-1]

    molecules_sorted_by_max_value = molecules_unsorted[sorted_indices]

    for i, molecule in enumerate(molecules_sorted_by_max_value):

        values = molecule_dict[molecule]

        if len(values) == 1:

            ax.plot(i, values[0], 'ko', ms=3)

        else:

            min_val = min(values)
            max_val = max(values)
            val_range = max_val-min_val

            ax.errorbar(i, (max_val+min_val)/2, yerr=val_range/2, fmt='None', ecolor='k', capsize=5)
    ax.set_xticks(np.arange(len(molecules_sorted_by_max_value)))
    ax.set_xticklabels([latex_molecule_name(mol) for mol in molecules_sorted_by_max_value], rotation=90)

    ax.set_ylabel(ylabel)

    return fig


def make_abundance_fraction_dict_from_coldens_range_dict(molecule_dict):

    # first we need to FIND the total amount of nitrogen, in coldens units...
    total_N_column = np.sum([np.max(y) for y in molecule_dict.values()])

    # then we need to divide all the values in the range dict by that number 
    # (non-destructively) and build a new dict of abundance fraction ranges.
    abundance_fraction_dict = {}

    for molecule in molecule_dict.keys():
        coldens_list = molecule_dict[molecule]
        abundance_fraction_dict[molecule] = [coldens/total_N_column for coldens in coldens_list]

    return abundance_fraction_dict


# included just for memory...
# contains_N = np.array(['N' in item for item in table['Molecule']])
# nitrogen_molecules_table = table[contains_N]

