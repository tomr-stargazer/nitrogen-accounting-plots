"""
Let's make a basic plot showing N/H2O values from protostars, comets, and ISM dust.

"""

import numpy as np
import matplotlib.pyplot as plt

# 0. Quick function defn.

def mean_and_yerr(min_val, max_val):
    mean = (min_val+max_val)/2
    yerr = np.abs(max_val-min_val)/2
    return mean, yerr

# 1. Numbers setup.

# Protostars:

okl_HCN = 0.10
e_okl_HCN = 0.06

okl_Norgs = 0.13
e_okl_Norgs = 0.07

iras_HCN = 1.2
e_iras_HCN = 0.8

# Cometary values:
min_HCN_comet = 0.07
max_HCN_comet = 0.45
HCN_comet, e_HCN_comet = mean_and_yerr(min_HCN_comet, max_HCN_comet)

min_NH3_comet = 0.3
max_NH3_comet = 1.5
NH3_comet, e_NH3_comet = mean_and_yerr(min_NH3_comet, max_NH3_comet)

upper_otherN_comet = 0.24

min_iceN_comet = 0.6
max_iceN_comet = 2.1
iceN_comet, e_iceN_comet = mean_and_yerr(min_iceN_comet, max_iceN_comet)

min_totN_comet = 6
max_totN_comet = 21
totN_comet, e_totN_comet = mean_and_yerr(min_totN_comet, max_totN_comet)

# ISM dust:
N_dust = 19
e_N_dust = 12


# 2. Plot setup

fig = plt.figure()
ax = fig.add_subplot(111)

x_protostar_offset = 0
x_orion = 1 +x_protostar_offset
x_iras = 2 +x_protostar_offset

x_comet_offset = 3 + x_protostar_offset
x_comet_HCN = 1 + x_comet_offset
x_comet_NH3 = 2 + x_comet_offset
x_comet_other = 3 + x_comet_offset
x_comet_ices = 2.1 + x_comet_offset
x_comet_total = 4 + x_comet_offset

x_ismdust_offset = 5 + x_comet_offset + x_protostar_offset
x_ismdust = 1 + x_ismdust_offset

# 3. Plot

# Orion
ax.errorbar(x_orion, okl_HCN, yerr=e_okl_HCN)
ax.errorbar(x_orion+0.1, okl_Norgs, yerr=e_okl_Norgs)

# IRAS
ax.errorbar(x_iras, iras_HCN, yerr=e_iras_HCN)

# Comet
ax.errorbar(x_comet_HCN, HCN_comet, yerr=e_HCN_comet)
ax.errorbar(x_comet_NH3, NH3_comet, yerr=e_NH3_comet)
ax.errorbar(x_comet_ices, iceN_comet, yerr=e_iceN_comet)
ax.errorbar(x_comet_total, totN_comet, yerr=e_totN_comet)
ax.plot(x_comet_other, upper_otherN_comet, 'v', ms=10)

# ISM Dust
ax.errorbar(x_ismdust, N_dust, yerr=e_N_dust)

# Vertical lines
ax.axvline(x_comet_offset, color='k', linestyle='--', lw=0.5)
ax.axvline(x_ismdust_offset, color='k', linestyle='--', lw=0.5)


plt.semilogy()
plt.show()
