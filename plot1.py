#!/usr/local/miniconda/bin/python
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
en = np.fromfile('en.dat')
flux = np.fromfile('flux.dat')
lum = flux * en * en * 1.65e42 * 1.60e-9

ax.set_ylabel(r'$\nu L_{\nu}\ {\rm[erg\ s^{-1}]}$')
ax.minorticks_on()
ax.set_xlabel('Energy (keV)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(5e35, 2e44)

plt.savefig('flux.pdf', bbox_inches='tight')
