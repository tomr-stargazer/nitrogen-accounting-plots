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


def coldens_range_dict_from_table(table):
    """ result: (main isotopologue -> tuple of coldens values) dict """

    # I need an intermediate table.
    # table with rows for each component -> table with rows for each species.
    adjusted_table = table.copy()
    adjusted_table.remove_columns(("theta_s", "T_rot", "v_lsr", "Deltav"))
    molecule_list = np.array(list(set(table['Molecule'])))
    molecule_count = np.array([np.sum(table['Molecule']==molecule) for molecule in molecule_list])
    duplicate_species = molecule_list[molecule_count>1]

    for molecule in duplicate_species:
        # get the 
        combined_N_tot = np.sum(table['N_tot'][table['Molecule'] == molecule])
        # do an operation wherein we remove the old ones
        indices = adjusted_table['Molecule'] == molecule
        adjusted_table.remove_rows(indices)
        adjusted_table.add_row((molecule, combined_N_tot))

    ratio_12C_13C = 69
    ratio_14N_15N = 388

    molecule_dict = {}

    for row in adjusted_table:
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


def abundance_fraction_dict_from_coldens_range_dict(molecule_dict):
    """ 
    This 'naive' approach is not quite right for the range of HCN values. 

    """

    # first we need to FIND the total amount of nitrogen, in coldens units...
    total_N_column = np.sum([np.median(y) for y in molecule_dict.values()])

    # then we need to divide all the values in the range dict by that number 
    # (non-destructively) and build a new dict of abundance fraction ranges.
    abundance_fraction_dict = {}

    for molecule in molecule_dict.keys():
        coldens_list = molecule_dict[molecule]
        abundance_fraction_dict[molecule] = [coldens/total_N_column for coldens in coldens_list]

    return abundance_fraction_dict


def abundance_fraction_dict_v2(molecule_dict):
    """
    This does a "poor-man's" error propagation to get the uncertainties under control.

    """

    naive_total_N_column = np.sum([np.median(y) for y in molecule_dict.values()])

    abundance_fraction_dict = {}

    for molecule in molecule_dict.keys():

        # This approach ensures that we never get an abundance > 1.
        N_tot_others = naive_total_N_column - np.median(molecule_dict[molecule])

        abundance_list = [value / (N_tot_others+value) for value in molecule_dict[molecule]]

        abundance_fraction_dict[molecule] = abundance_list

    return abundance_fraction_dict


def make_charts_for_OrionKL_HotCore():

    table = astropy.table.Table.read("apj493753t4_ascii_copy_hotcore.txt", format='ascii.basic', delimiter='\t', guess=False, data_start=3)    
    contains_N = np.array(['N' in item for item in table['Molecule']])
    nitrogen_molecules_table = table[contains_N]
    molecule_dict = coldens_range_dict_from_table(nitrogen_molecules_table)
    del molecule_dict['NH_2D']
    del molecule_dict['DCN']
    fraction_dict_v2 = abundance_fraction_dict_v2(molecule_dict)

    figure_list = []
    figure_list.append(make_errorbar_figure(molecule_dict))
    plt.semilogy()
    plt.savefig("column_density_plot.png", bbox_inches='tight')

    figure_list.append(make_errorbar_figure(fraction_dict_v2, ylabel="Nitrogen fraction"))
    plt.semilogy()
    plt.savefig("nitrogen_fraction_plot.png", bbox_inches='tight')

    hydrogen_column_density = 3.1e23
    X_abundance_dict = {x: np.array(y)/hydrogen_column_density for x, y in molecule_dict.items()}

    figure_list.append(make_errorbar_figure(X_abundance_dict, ylabel="Molecular abundance X"))
    plt.semilogy()
    plt.savefig("abundance_plot.png", bbox_inches='tight')

    Si_to_H_ratio = 3.16e-05
    N_to_Si_conversion_constant = 1/2 * 1/Si_to_H_ratio
    N_Si_ratio_dict = {x: y*N_to_Si_conversion_constant for x, y in X_abundance_dict.items()}

    figure_list.append(make_errorbar_figure(N_Si_ratio_dict, ylabel="N/Si ratio for each molecule"))
    plt.semilogy()
    plt.savefig("N_Si_plot.png", bbox_inches='tight')

    return figure_list

