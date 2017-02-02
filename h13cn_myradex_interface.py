"""
This is a boilerplate-skipping wrapper for myradex h13cn.

It's not very gerneralizable.

"""

from __future__ import division
import wrapper_my_radex
wrapper = wrapper_my_radex.myradex_wrapper
from myradex_table_function import table_from_myradex_datatransitions

n_levels, n_item, n_transitions = wrapper.config_basic(
    '/Users/tsrice/Documents/Academia/Misc_Software/Radex/data/',
    'h13cn@xpol.dat', 2.7, True)

default_params = {
    'tkin': None,
    'dv_cgs': 1e5, # What is this? oh, it's delta v in cm/scp
    'dens_x_cgs': 1e-2, # I think this is not a relevant parameter?
    'ncol_x_cgs': None,
    'h2_density_cgs': None,
    'hi_density_cgs': 0.0,
    'oh2_density_cgs': 0.0,
    'ph2_density_cgs': 0.0,
    'hii_density_cgs': 0.0,
    'electron_density_cgs': 0.0,
    'n_levels': n_levels,
    'n_item': n_item,
    'n_transitions': n_transitions,
    'geotype': 'lvg'}

# How do we pass in that dict of default params but update specific items if necessary?
def myradex_h13cn(kinetic_temperature=None, collider_density=None, column_density=None):

    params = default_params.copy()

    params['tkin'] = kinetic_temperature
    params['ncol_x_cgs'] = column_density
    params['h2_density_cgs'] = collider_density

    # load the stuff from the molecular data file

    output_dict = {}
    (output_dict['energies'], 
        output_dict['f_occupations'], 
        output_dict['data_transitions'], 
        output_dict['cooling_rate']) = wrapper.run_one_params(**params)
    output_dict['table'] = table_from_myradex_datatransitions(output_dict['data_transitions'])

    return output_dict
