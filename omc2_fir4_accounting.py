"""
Here I'm gonna do some stuff.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as c

# First: the numbers.

# Crudest possible thing... let's just dump all the X numbers into one array.

omc2_fir4_X_HCN_values = np.array([
	9.7e-9, 44e-9, 2.3e-9, 11e-9, 5.2e-9, 23e-9, 1.2e-9, 5.8e-9])

omc2_fir4_X_HCN_from_H13CN = np.array([
	9.7e-9, 44e-9, 5.2e-9, 23e-9])

omc2_fir4_X_HCN_from_HC15N = np.array([
	2.3e-9, 11e-9, 1.2e-9, 5.8e-9])


fig = plt.figure()

plt.plot(np.ones_like(omc2_fir4_X_HCN_from_H13CN)+0.1, omc2_fir4_X_HCN_from_H13CN, 'ro')
plt.plot(np.ones_like(omc2_fir4_X_HCN_from_HC15N)-0.1, omc2_fir4_X_HCN_from_HC15N, 'bo')

plt.ylabel("$X_{HCN}$ derived from HCN isotopologues")

plt.show()

# Now let's divide these by the H2/Si ratio... digging it up from Ted's stuff now.

nH_per_H2 = 2

# The cosmic abundance of Si per H nucleus is... (see Nieva & Przybilla 2012)
nSi_per_H = 28e-6

omc2_fir4_N_Si_values_13C = omc2_fir4_X_HCN_from_H13CN / nH_per_H2 / nSi_per_H
omc2_fir4_N_Si_values_15N = omc2_fir4_X_HCN_from_HC15N / nH_per_H2 / nSi_per_H


omc2_fir4_N_Si_values_13C = omc2_fir4_X_HCN_from_H13CN / nH_per_H2 / nSi_per_H
omc2_fir4_N_Si_values_15N = omc2_fir4_X_HCN_from_HC15N / nH_per_H2 / nSi_per_H

fig2 = plt.figure()

plt.plot(np.ones_like(omc2_fir4_X_HCN_from_H13CN)+0.1, omc2_fir4_N_Si_values_13C, 'ro-')
plt.plot(np.ones_like(omc2_fir4_X_HCN_from_HC15N)-0.1, omc2_fir4_N_Si_values_15N, 'bo-')

# Values derived from Orion KL analysis:
N_Si_value_from_KL_min = 0.02
N_Si_value_from_KL_max = 0.028

N_Si_IRAS2a = 1e-3

plt.plot([2, 2], [N_Si_value_from_KL_min, N_Si_value_from_KL_max], 'g', lw=3)
plt.plot(3, N_Si_IRAS2a, 'ko')

plt.ylabel("$n_{HCN}/n_{Si}$", fontsize=18)
plt.NullFormatter())

plt.semilogy()
plt.xlim(0.5, 3.5)
plt.title("Comparison of N/Si between protostar regions")

plt.text(1.1, 1e-2, "HCN from H$^{13}$CN", color='r', rotation='vertical')
plt.text(0.9, 3e-3, "HCN from HC$^{15}$N", color='b', rotation='vertical')
plt.text(0.7, 3e-3, "OMC-2 FIR4 (Shimajiri+15)", rotation='vertical')
plt.text(2.9, 0.001, "NGC 1333 IRAS 2A", rotation='vertical', verticalalignment='center')
plt.text(1.9, 0.02, "Orion KL", rotation='vertical', verticalalignment='center')

plt.show()


