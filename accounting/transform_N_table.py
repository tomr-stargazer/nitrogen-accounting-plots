"""
Doing some quick transformations of the nitrogen table.

"""

from __future__ import division
import numpy as np
import astropy.table

def transform_table(table, species=('hcn', 'nh3')):
    """ 
    Format of table:
    Molecule theta_s T_rot N_tot v_lsr Deltav

    first col is strings, all others are floats

    """

    new_table = table.copy()

    # 1. the HCN we want to show is the HC15N value, scaled up (via isotopic ratios) to HCN. 

    if 'hcn' in species:
        # so -
        # b. hide the other HCN values
        # remove some rows
        contains_hcn_fn = lambda x: ('HCN' in x) or ('DCN' in x) or ('H^13CN' in x) or ('HC^15N' in x)
        hcn_row_indices = np.array([contains_hcn_fn(row['Molecule']) for row in new_table])
        new_table.remove_rows(hcn_row_indices)
        # a. make the "real" value we want to show
        # make it
        ratio_12C_13C = 69 # from Wilson '99 local ISM
        real_hcn_N_tot = np.sum(table['N_tot'][table['Molecule'] == 'H^13CN,nu_2 = 1']) * ratio_12C_13C
        real_hcn_row = ('HCN', 0, 0, real_hcn_N_tot, 0, 0)
        # add it as a row
        new_table.add_row(real_hcn_row)

    if 'nh3' in species:    
        # 2. Doing the same thing with NH3.
        # b. hide the other HCN values
        # remove some rows
        contains_nh3_fn = lambda x: ('NH_3' in x) or ('^15NH_3' in x) or ('NH_2D' in x)
        nh3_row_indices = np.array([contains_nh3_fn(row['Molecule']) for row in new_table])
        new_table.remove_rows(nh3_row_indices)
        # a. make the "real" value we want to show
        # make it
        ratio_14n_15n = 388
        real_nh3_N_tot = np.sum(table['N_tot'][table['Molecule'] == '^15NH_3']) * ratio_14n_15n
        real_nh3_row = ('NH_3', 0, 0, real_nh3_N_tot, 0, 0)
        # add it as a row
        new_table.add_row(real_nh3_row)

    return new_table




