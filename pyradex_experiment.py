""" I am experimenting with pyradex. """

from __future__ import division

import numpy as np
import matplotlib as plt
import pyradex

# how do I work with pyradex?

R = pyradex.Radex(species='c18o', abundance=1e-8, column=1e14, temperature=20)
R.run_radex()

# Ok CAN I DO A THING

T = R(collider_densities={'H2':1000}, column=10**12)

T.pprint()

# tom, let's do the shameful thing and do a LOOP

rho_array = np.logspace(1, 6.5, 20)
Tex_array_1_0 = np.zeros_like(rho_array)
Tex_array_2_1 = np.zeros_like(rho_array)
Tex_array_3_2 = np.zeros_like(rho_array)

for i, rho in enumerate(np.logspace(1, 6.5, 20)):

    new_T = R(collider_densities={'H2':rho})

    Tex_array_1_0[i] = new_T[0]['Tex']
    Tex_array_2_1[i] = new_T[1]['Tex']
    Tex_array_3_2[i] = new_T[2]['Tex']

print Tex_array_1_0
print Tex_array_2_1
print Tex_array_3_2

fig = plt.Figure()

plt.plot(np.log10(rho_array), Tex_array_1_0)
plt.plot(np.log10(rho_array), Tex_array_2_1)
plt.plot(np.log10(rho_array), Tex_array_3_2)

plt.ylim(0,25)
plt.xlim(1,6.5)

plt.show()


