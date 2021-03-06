"""
Making plots for h13cn (and friends), using the
goldsmith/plume plots as a model.

"""

from __future__ import division
import time
import pdb

import numpy as np
import matplotlib.pyplot as plt
import pyradex

import h13cn_myradex_interface

# This import will have to change if I refactor correction_factor much.
from correction_factor.correction_factor_2 import correction_factor_given_radex_table, correction_factor_given_myradex_output

# h13cn_J_upper_list = [1, 3, 4, 6, 7, 8, 9, 10]
h13cn_J_upper_list = [3, 4, 6, 7, 8, 9, 10]
h13cn_J_lower_list = [x-1 for x in h13cn_J_upper_list]

def h13cn_Tex_vs_density_A(save=True, column=1e16, print_timing=False):
    """ 
    Generates a clone of Fig. 9 from Goldsmith+'97 for h13cn.

    Parameters chosen for similarity with Goldsmith '97.

    """

    start = time.time()

    # Here we have to specify 2 of 3: abundance, temperature, density.
    R = pyradex.Radex(species='h13cn@xpol', density=1e1, temperature=20, column=column)

    fig = plt.figure()

    ax = fig.add_subplot(111)

    density_array = np.logspace(1, 6.5, 100)

    Tex_array_1_0 = np.zeros_like(density_array)
    Tex_array_3_2 = np.zeros_like(density_array)
    Tex_array_4_3 = np.zeros_like(density_array)

    for i, density in enumerate(density_array):

        new_T = R(collider_densities={'H2':density})

        Tex_array_1_0[i] = new_T[0]['Tex']
        Tex_array_3_2[i] = new_T[2]['Tex']
        Tex_array_4_3[i] = new_T[3]['Tex']

    ax.plot(np.log10(density_array), Tex_array_1_0, '-', lw=1.5, color='b')
    ax.plot(np.log10(density_array), Tex_array_3_2, '--', lw=1.5, color='r')
    ax.plot(np.log10(density_array), Tex_array_4_3, '-.', lw=2, color='purple')

    ax.set_ylim(0,25)
    ax.set_xlim(1,6.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("EXCITATION TEMPERATURE/K", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(1.5, 20, "H$^{13}$CN, T = 20 K", fontsize=18, family='serif')
    ax.text(5, 22, "$J = 1 - 0$", fontsize=16, family='serif', color='b')
    ax.text(3.6, 7, "$J = 4 - 3$", fontsize=16, family='serif', color='purple')
    ax.text(5.6, 4, "$J = 3 - 2$", fontsize=16, family='serif', color='r')

    ax.set_title("version A")
    ax.minorticks_on()

    if save:
        fig.savefig("h13cn_Tex_vs_density_A.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig


def h13cn_Tex_vs_density_A_myradex(save=True, column=1e16, print_timing=False):
    """ 
    Generates a clone of Fig. 9 from Goldsmith+'97 for h13cn.

    Parameters chosen for similarity with Goldsmith '97.

    """

    start = time.time()

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_array = np.logspace(1, 6.5, 100)

    Tex_array_1_0 = np.zeros_like(density_array)
    Tex_array_3_2 = np.zeros_like(density_array)
    Tex_array_4_3 = np.zeros_like(density_array)

    for i, density in enumerate(density_array):

        od = h13cn_myradex_interface.myradex_h13cn(kinetic_temperature=20, column_density=column, collider_density=density)
        new_T = od['table']

        Tex_array_1_0[i] = new_T[0]['Tex']
        Tex_array_3_2[i] = new_T[2]['Tex']
        Tex_array_4_3[i] = new_T[3]['Tex']

    ax.plot(np.log10(density_array), Tex_array_1_0, '-', lw=1.5, color='b')
    ax.plot(np.log10(density_array), Tex_array_3_2, '--', lw=1.5, color='r')
    ax.plot(np.log10(density_array), Tex_array_4_3, '-.', lw=2, color='purple')

    ax.set_ylim(0,25)
    ax.set_xlim(1,6.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("EXCITATION TEMPERATURE/K", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(1.5, 20, "H$^{13}$CN, T = 20 K", fontsize=18, family='serif')
    ax.text(5, 22, "$J = 1 - 0$", fontsize=16, family='serif', color='b')
    ax.text(3.6, 7, "$J = 4 - 3$", fontsize=16, family='serif', color='purple')
    ax.text(5.6, 4, "$J = 3 - 2$", fontsize=16, family='serif', color='r')

    ax.set_title("version A.myradex")
    ax.minorticks_on()

    if save:
        fig.savefig("h13cn_Tex_vs_density_A_myradex.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig


def h13cn_Tex_vs_density_B(save=True, column=1e16, print_timing=False):
    """ 
    Generates a clone of Fig. 9 from Goldsmith+'97 for h13cn.

    Parameters chosen for relevance to astrophysical hot corino.

    """

    start = time.time()

    # Here we have to specify 2 of 3: abundance, temperature, density.
    R = pyradex.Radex(species='h13cn@xpol', density=1e4, temperature=100, column=column)

    fig = plt.figure()

    ax = fig.add_subplot(111)

    density_array = np.logspace(4, 8.5, 100)

    Tex_array_1_0 = np.zeros_like(density_array)
    Tex_array_3_2 = np.zeros_like(density_array)
    Tex_array_4_3 = np.zeros_like(density_array)

    for i, density in enumerate(density_array):

        new_T = R(collider_densities={'H2':density})

        Tex_array_1_0[i] = new_T[0]['Tex']
        Tex_array_3_2[i] = new_T[2]['Tex']
        Tex_array_4_3[i] = new_T[3]['Tex']

    ax.plot(np.log10(density_array), Tex_array_1_0, '-', lw=1.5, color='b')
    ax.plot(np.log10(density_array), Tex_array_3_2, '--', lw=1.5, color='r')
    ax.plot(np.log10(density_array), Tex_array_4_3, '-.', lw=2, color='purple')

    ax.set_ylim(0,130)
    ax.set_xlim(4,8.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("EXCITATION TEMPERATURE/K", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(4.1, 100, "H$^{13}$CN, T = 100 K", fontsize=18, family='serif')
    ax.text(5, 80, "$J = 1 - 0$", fontsize=16, family='serif', color='b')
    ax.text(5, 15, "$J = 4 - 3$", fontsize=16, family='serif', color='purple')
    ax.text(6, 80, "$J = 3 - 2$", fontsize=16, family='serif', color='r')

    ax.set_title("version B")
    ax.minorticks_on()

    if save:
        fig.savefig("h13cn_Tex_vs_density_B.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig


def h13cn_fc_vs_density_A(save=True, print_timing=False):
    """ 
    Generates a clone of Fig. 10 from Goldsmith+'97 for h13cn 

    Parameters chosen for similarity with Goldsmith '97.

    """

    start = time.time()

    # Here we have to specify 2 of 3: abundance, temperature, density.
    R = pyradex.Radex(species='h13cn@xpol', density=1e1, temperature=20, column=1e16)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_array = np.logspace(1, 7, 100)
    temperature_list = [10, 30, 50]

    f_c_array_10K = np.zeros_like(density_array)
    f_c_array_30K = np.zeros_like(density_array)
    f_c_array_50K = np.zeros_like(density_array)

    f_c_array_list = [f_c_array_10K, f_c_array_30K, f_c_array_50K]

    for j, temp in enumerate(temperature_list):

        R.temperature = temp

        for i, rho in enumerate(density_array):

            new_T = R(collider_densities={'H2':rho})

            # use the correction_factor code here!
            f_c = correction_factor_given_radex_table(new_T, h13cn_J_lower_list)

            f_c_array_list[j][i] = f_c
            
    ax.plot(np.log10(density_array), f_c_array_10K, '-', lw=1.5)
    ax.plot(np.log10(density_array), f_c_array_30K, ':', lw=2)
    ax.plot(np.log10(density_array), f_c_array_50K, '--', lw=1.5)

    ax.set_ylim(1, 3)
    ax.set_xlim(2, 7.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("CORRECTION FACTOR FOR COLUMN DENSITY", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(3, 2, "H$^{13}$CN", fontsize=18, family='serif')
    ax.text(4.5, 2.5, r"$T_k = 10\ \rm{K}$", fontsize=16, family='serif', color='b')
    ax.text(4.6, 2.25, r"$T_k = 30\ \rm{K}$", fontsize=16, family='serif', color='g')
    ax.text(6, 1.5, r"$T_k = 50\ \rm{K}$", fontsize=16, family='serif', color='r')
    ax.text(2.5, 1.2, r"$J_u=1,\ 3,\ 4,\ 6,\ 7,\ 8,\ 9,\ 10$ OBSERVED", fontsize=14, family='serif')

    ax.set_title("version A")
    ax.minorticks_on()

    if save:
        fig.savefig("h13cn_fc_vs_density_A.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig


def h13cn_fc_vs_density_A_myradex(save=True, print_timing=False):
    """ 
    Generates a clone of Fig. 10 from Goldsmith+'97 for h13cn 

    Parameters chosen for similarity with Goldsmith '97.

    """

    start = time.time()

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_array = np.logspace(1, 7, 100)
    temperature_list = [10, 30, 50]
    column=1e10

    f_c_array_10K = np.zeros_like(density_array)
    f_c_array_30K = np.zeros_like(density_array)
    f_c_array_50K = np.zeros_like(density_array)

    f_c_array_list = [f_c_array_10K, f_c_array_30K, f_c_array_50K]

    for j, temp in enumerate(temperature_list):

        for i, rho in enumerate(density_array):

            od = h13cn_myradex_interface.myradex_h13cn(kinetic_temperature=temp, 
                column_density=column, collider_density=rho)
            new_T = od['table']

            # use the correction_factor code here!
            f_c = correction_factor_given_myradex_output(od['f_occupations'], h13cn_J_lower_list)

            f_c_array_list[j][i] = f_c
            
    ax.plot(np.log10(density_array), f_c_array_10K, '-', lw=1.5)
    ax.plot(np.log10(density_array), f_c_array_30K, ':', lw=2)
    ax.plot(np.log10(density_array), f_c_array_50K, '--', lw=1.5)

    ax.set_ylim(1, 3)
    ax.set_xlim(2, 7.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("CORRECTION FACTOR FOR COLUMN DENSITY", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(3, 2, "H$^{13}$CN", fontsize=18, family='serif')
    ax.text(4.5, 2.5, r"$T_k = 10\ \rm{K}$", fontsize=16, family='serif', color='b')
    ax.text(4.6, 2.25, r"$T_k = 30\ \rm{K}$", fontsize=16, family='serif', color='g')
    ax.text(6, 1.5, r"$T_k = 50\ \rm{K}$", fontsize=16, family='serif', color='r')
    ax.text(2.5, 1.2, r"$J_u=1,\ 3,\ 4,\ 6,\ 7,\ 8,\ 9,\ 10$ OBSERVED", fontsize=14, family='serif')

    ax.set_title("version A.myradex")
    ax.minorticks_on()

    if save:
        fig.savefig("h13cn_fc_vs_density_A_myradex.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig


def h13cn_fc_vs_density_B(save=True, print_timing=False):
    """ 
    Generates a clone of Fig. 10 from Goldsmith+'97 for h13cn 

    Parameters chosen for similarity with astrophysical hot corino.

    """

    start = time.time()

    # Here we have to specify 2 of 3: temperature, density.
    R = pyradex.Radex(species='h13cn@xpol', density=1e4, temperature=100, column=1e10, escapeProbGeom='lvg')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_array = np.logspace(4, 8.5, 100)
    temperature_list = [50, 100, 200]

    f_c_array_50K = np.zeros_like(density_array)
    f_c_array_100K = np.zeros_like(density_array)
    f_c_array_200K = np.zeros_like(density_array)

    f_c_array_list = [f_c_array_50K, f_c_array_100K, f_c_array_200K]

    for j, temp in enumerate(temperature_list):

        R.temperature = temp

        for i, rho in enumerate(density_array):

            new_T = R(collider_densities={'H2':rho})

            # use the correction_factor code here!
            f_c = correction_factor_given_radex_table(new_T, h13cn_J_lower_list)

            f_c_array_list[j][i] = f_c
            
    ax.plot(np.log10(density_array), f_c_array_50K, 'b-', lw=1.5)
    ax.plot(np.log10(density_array), f_c_array_100K, 'g:', lw=2)
    ax.plot(np.log10(density_array), f_c_array_200K, 'r--', lw=1.5)

    ax.set_ylim(1, 2.5)
    ax.set_xlim(4, 8.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("CORRECTION FACTOR FOR COLUMN DENSITY", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(4.5, 2, "H$^{13}$CN", fontsize=18, family='serif')
    ax.text(6.75, 1.65, r"$T_k = 50\ \rm{K}$", fontsize=16, family='serif', color='b')
    ax.text(7.5, 1.3, r"$T_k = 100\ \rm{K}$", fontsize=16, family='serif', color='g')
    ax.text(7.75, 1.75, r"$T_k = 200\ \rm{K}$", fontsize=16, family='serif', color='r')
    ax.text(4.5, 1.2, r"$J_u=1,\ 3,\ 4,\ 6,\ 7,\ 8,\ 9,\ 10$ OBSERVED", fontsize=14, family='serif')

    ax.set_title("version B")
    ax.minorticks_on()

    if save:
        fig.savefig("h13cn_fc_vs_density_B.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig


def h13cn_fc_vs_density_B_myradex(save=True, print_timing=False):
    """ 
    Generates a clone of Fig. 10 from Goldsmith+'97 for h13cn 

    Parameters chosen for similarity with astrophysical hot corino.

    """

    start = time.time()

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_array = np.logspace(4, 8.5, 100)
    temperature_list = [50, 100, 200]
    column=1e10

    f_c_array_50K = np.zeros_like(density_array)
    f_c_array_100K = np.zeros_like(density_array)
    f_c_array_200K = np.zeros_like(density_array)

    f_c_array_list = [f_c_array_50K, f_c_array_100K, f_c_array_200K]

    for j, temp in enumerate(temperature_list):

        for i, rho in enumerate(density_array):

            od = h13cn_myradex_interface.myradex_h13cn(kinetic_temperature=temp, 
                column_density=column, collider_density=rho)
            new_T = od['table']

            # use the correction_factor code here!
            f_c = correction_factor_given_myradex_output(od['f_occupations'], h13cn_J_lower_list)

            f_c_array_list[j][i] = f_c
            
    ax.plot(np.log10(density_array), f_c_array_50K, 'b-', lw=1.5)
    ax.plot(np.log10(density_array), f_c_array_100K, 'g:', lw=2)
    ax.plot(np.log10(density_array), f_c_array_200K, 'r--', lw=1.5)

    ax.set_ylim(1, 2.5)
    ax.set_xlim(4, 8.5)

    fontdict = {'family':'serif'}

    ax.set_ylabel("CORRECTION FACTOR FOR COLUMN DENSITY", fontdict=fontdict)
    ax.set_xlabel("LOG (MOLECULAR HYDROGEN DENSITY/CM$^{-3}$)", fontdict=fontdict)

    ax.text(4.5, 2, "H$^{13}$CN", fontsize=18, family='serif')
    ax.text(6.75, 1.65, r"$T_k = 50\ \rm{K}$", fontsize=16, family='serif', color='b')
    ax.text(7.5, 1.3, r"$T_k = 100\ \rm{K}$", fontsize=16, family='serif', color='g')
    ax.text(7.75, 1.75, r"$T_k = 200\ \rm{K}$", fontsize=16, family='serif', color='r')
    string_of_J_uppers = r'\ ,'.join([str(x) for x in h13cn_J_upper_list])
    ax.text(4.5, 1.2, r"$J_u={0}$ OBSERVED".format(string_of_J_uppers), fontsize=14, family='serif')

    ax.set_title("version B.myradex")
    ax.minorticks_on()

    if save:
        fig.savefig("h13cn_fc_vs_density_B_myradex.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig


def h13cn_fc_vs_temperature_A(save=True, print_timing=False, n_points=20, column=5e16, deltav=1):
    """ 
    Generates a clone of Fig. 3 from Plume+'12 for h13cn 

    Parameters chosen to match Plume '12.

    """

    # Plume '12 claims to use an abundance of 2e-7 (i.e. one C18O for every 500 x 10^4 H2 atoms)
    # but I find that 5e-7 matches the appearance much more closely.

    start = time.time()

    # J_upper_list = [2,3,5,7,8,9,11,15]
    # J_lower_list = [x-1 for x in J_upper_list]

    # Here we have to specify 2 of 3: abundance, temperature, density.
    # R = pyradex.Radex(species='h13cn@xpol', abundance=abundance, column=column, temperature=10, escapeProbGeom='lvg', deltav=deltav)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_list = [1e5, 1e6, 1e7]
    temperature_array = np.linspace(10, 400, n_points)

    f_c_array_n5 = np.zeros_like(temperature_array)
    f_c_array_n6 = np.zeros_like(temperature_array)
    f_c_array_n7 = np.zeros_like(temperature_array)

    f_c_array_list = [f_c_array_n5, f_c_array_n6, f_c_array_n7]

    # for j, density in enumerate(density_list):

    #     R.density = density

    for i, temp in enumerate(temperature_array):
        # R.temperature = temp

        for j, density in enumerate(density_list):

            # R.density = density
            # R.temperature = temp
            # R.column = column

            try:
                del R
            except:
                pass
            R = pyradex.Radex(species='h13cn@xpol', column=column, density=density, temperature=temp, escapeProbGeom='sphere')
            R.run_radex(reuse_last=False, reload_molfile=True)

            new_T = R.get_table()

            # use the correction_factor code here!
            f_c = correction_factor_given_radex_table(new_T, h13cn_J_lower_list)
            # They only use the lowest 18 rotational levels, 
            # even though the plot's appearance changes when 40 are included.

            if np.sum(new_T['lowerlevelpop']) < 0.99:
                print( "unphysical population statistics encountered")
                print( "conditions: T={0}, n={1}".format(temp, density))
                print( "abundance: {0}, coldens: {1}, deltav={2}".format(abundance, column, deltav))
                f_c_array_list[j][i] = np.nan
                print(r"run %debug and explore the variable `new_T` here.")
                pdb.set_trace()
            elif min(new_T['Tex']) < 0:
                print( "negative Tex encountered")
                print( "conditions: T={0}, n={1}".format(temp, density))
                print( "abundance: {0}, coldens: {1}, deltav={2}".format(abundance, column, deltav))
                # raise ValueError("intentionally erroring")
                f_c_array_list[j][i] = f_c
                pdb.set_trace()
            else:
                f_c_array_list[j][i] = f_c
            
    ax.plot(temperature_array, f_c_array_n7, 'b-', lw=1.5, label=r'$n = 10^7\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n6, 'g--', lw=1.5, label=r'$n = 10^6\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n5, 'r:', lw=2, label=r'$n = 10^5\ \rm{cm}^{-3}$')

    ax.legend(frameon=False, loc='upper center', fontsize=18)

    ax.set_xlim(0, 400)
    ax.set_ylim(1.3, 2)
    # ax.set_xticks([0,100,200,300,400])
    # ax.set_yticks([1.6,1.7,1.8,1.9,2])

    fontdict = {'family':'serif', 'fontsize':16}

    ax.set_ylabel("Correction Factor", fontdict=fontdict)
    ax.set_xlabel("Temp (K)", fontdict=fontdict)

    ax.minorticks_on()

    if save:
        fig.savefig("h13cn_fc_vs_temperature_A.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig    


def h13cn_fc_vs_temperature_B(save=True, print_timing=False, n_points=20):
    """ 
    Generates a clone of Fig. 3 from Plume+'12 for h13cn 

    Parameters chosen to match astrophysical hot corino.

    """

    # Plume '12 claims to use an abundance of 2e-7 (i.e. one C18O for every 500 x 10^4 H2 atoms)
    # but I find that 5e-7 matches the appearance much more closely.

    start = time.time()

    # J_upper_list = [2,3,5,7,8,9,11,15]
    # J_lower_list = [x-1 for x in J_upper_list]

    # Here we have to specify 2 of 3: abundance, temperature, density.
    R = pyradex.Radex(species='h13cn@xpol', density=1e6, column=1e10, temperature=50, escapeProbGeom='lvg')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_list = [1e6, 1e7, 1e8, 1e9]
    temperature_array = np.linspace(10, 400, n_points)

    f_c_array_n6 = np.zeros_like(temperature_array)
    f_c_array_n7 = np.zeros_like(temperature_array)
    f_c_array_n8 = np.zeros_like(temperature_array)
    f_c_array_n9 = np.zeros_like(temperature_array)

    f_c_array_list = [f_c_array_n6, f_c_array_n7, f_c_array_n8, f_c_array_n9]

    for j, density in enumerate(density_list):

        R.density = density

        for i, temp in enumerate(temperature_array):

            new_T = R(temperature=temp)

            # use the correction_factor code here!
            f_c = correction_factor_given_radex_table(new_T, h13cn_J_lower_list)

            f_c_array_list[j][i] = f_c
            
    ax.plot(temperature_array, f_c_array_n6, 'g--', lw=1.5, label=r'$n = 10^6\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n7, 'b-', lw=1.5, label=r'$n = 10^7\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n8, 'r:', lw=2, label=r'$n = 10^8\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n9, '-.', color='purple', lw=2, label=r'$n = 10^9\ \rm{cm}^{-3}$')

    ax.legend(frameon=False, loc='upper center', fontsize=18)

    ax.set_xlim(0, 400)
    # ax.set_ylim(1.6, 2)
    # ax.set_xticks([0,100,200,300,400])
    # ax.set_yticks([1.6,1.7,1.8,1.9,2])

    fontdict = {'family':'serif', 'fontsize':16}

    ax.set_ylabel("Correction Factor", fontdict=fontdict)
    ax.set_xlabel("Temp (K)", fontdict=fontdict)

    ax.minorticks_on()

    if save:
        fig.savefig("h13cn_fc_vs_temperature_B.pdf")

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig        

def h13cn_fc_vs_temperature_myradex(save=True, print_timing=False, n_points=20):
    """ 
    Generates a clone of Fig. 3 from Plume+'12 for h13cn 

    Parameters chosen to match astrophysical hot corino.

    """

    # Plume '12 claims to use an abundance of 2e-7 (i.e. one C18O for every 500 x 10^4 H2 atoms)
    # but I find that 5e-7 matches the appearance much more closely.

    column=1e10
    start = time.time()

    # J_upper_list = [2,3,5,7,8,9,11,15]
    # J_lower_list = [x-1 for x in J_upper_list]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    density_list = [1e6, 1e7, 1e8, 1e9]
    temperature_array = np.linspace(10, 400, n_points)

    f_c_array_n6 = np.zeros_like(temperature_array)
    f_c_array_n7 = np.zeros_like(temperature_array)
    f_c_array_n8 = np.zeros_like(temperature_array)
    f_c_array_n9 = np.zeros_like(temperature_array)

    f_c_array_list = [f_c_array_n6, f_c_array_n7, f_c_array_n8, f_c_array_n9]

    for j, density in enumerate(density_list):

        for i, temp in enumerate(temperature_array):

            od = h13cn_myradex_interface.myradex_h13cn(kinetic_temperature=temp, 
                column_density=column, collider_density=density)
            new_T = od['table']

            # use the correction_factor code here!
            f_c = correction_factor_given_myradex_output(od['f_occupations'], h13cn_J_lower_list)

            f_c_array_list[j][i] = f_c
            
    ax.plot(temperature_array, f_c_array_n6, 'g--', lw=1.5, label=r'$n = 10^6\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n7, 'b-', lw=1.5, label=r'$n = 10^7\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n8, 'r:', lw=2, label=r'$n = 10^8\ \rm{cm}^{-3}$')
    ax.plot(temperature_array, f_c_array_n9, '-.', color='purple', lw=2, label=r'$n = 10^9\ \rm{cm}^{-3}$')

    ax.legend(frameon=False, loc='upper center', fontsize=18)

    ax.set_xlim(0, 400)
    # ax.set_ylim(1.6, 2)
    # ax.set_xticks([0,100,200,300,400])
    # ax.set_yticks([1.6,1.7,1.8,1.9,2])

    fontdict = {'family':'serif', 'fontsize':16}

    ax.set_ylabel("Correction Factor", fontdict=fontdict)
    ax.set_xlabel("Temp (K)", fontdict=fontdict)

    ax.minorticks_on()

    if save:
        ax.set_ylim(1, 4)
        fig.savefig("h13cn_fc_vs_temperature_myradex.pdf")
        np.savetxt("h13cn_fc_v_T_myr_temperature_array.txt", temperature_array)
        np.savetxt("h13cn_fc_v_T_myr_fc_array_n6.txt", f_c_array_n6)
        np.savetxt("h13cn_fc_v_T_myr_fc_array_n7.txt", f_c_array_n7)
        np.savetxt("h13cn_fc_v_T_myr_fc_array_n8.txt", f_c_array_n8)
        np.savetxt("h13cn_fc_v_T_myr_fc_array_n9.txt", f_c_array_n9)

    plt.show()

    end = time.time()
    if print_timing:
        print (end-start)

    return fig            
