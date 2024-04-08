#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

from ase.io import read
from ase.calculators.aims import Aims
from ase.calculators.vasp import Vasp

from optparse import OptionParser
import os
import sys

usage = """%prog [options] xxx.cif"""
parser = OptionParser(usage)
(options, args) = parser.parse_args()
atoms= read(args[0], format='vasp')

encuts = [encut for encut in range(600, 1000, 50)] # 200~800
energies = []
for index, encut in enumerate(encuts):
    vasp_directory = '{:0>2d}.encut-{}'.format(index, encut)
    if not os.path.exists(vasp_directory):
        os.makedirs(vasp_directory)
    os.chdir(vasp_directory)
    calculator = Vasp(xc='PBE', kpts=(6, 6, 6),
                    encut=encut,
                    ismear=-5,
                    #luse_vdw=True,
                    lreal=False, lcharg=False, lwave=False
                    )
    atoms.set_calculator(calculator)
    energy = atoms.get_potential_energy()
    energies.append(energy)
    print('{:20s}: {:10.3f}'.format(vasp_directory, energy))
    os.chdir(os.path.dirname(os.getcwd())) # back to root directory
plt.plot(encuts, energies, '-o')
plt.xlabel('Cutoff Energy (eV)')
plt.ylabel('Total Energy (eV)')
plt.title('PBE encut test')
plt.savefig('encut_convergence.png')
plt.close()
