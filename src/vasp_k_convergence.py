import os
import numpy as np
import matplotlib.pyplot as plt
from ase.calculators.vasp import Vasp

from optparse import OptionParser
import os
import sys

usage = """%prog [options] xxx.cif"""
parser = OptionParser(usage)
(options, args) = parser.parse_args()
atoms= read(args[0], format='vasp')

k_points = [k for k in range(6, 20, 1)]
energies = []
for index, k in enumerate(k_points):
    vasp_directory = '{:0>2d}.kpt-{}'.format(index, k)
    if not os.path.exists(vasp_directory):
        os.makedirs(vasp_directory)
    os.chdir(vasp_directory)
    calculator = Vasp(xc='PBE', kpts=(k, k, k),
                    encut=600,
                    ismear=-5,
                    lreal=False, lcharg=False, lwave=False
                    )
    atoms.set_calculator(calculator)
    energy = atoms.get_potential_energy()
    energies.append(energy)
    print('{:20s}: {:10.3f}'.format(vasp_directory, energy))
    os.chdir(os.path.dirname(os.getcwd())) # back to root directory
plt.plot(k_points, energies, '-o')
plt.xlabel('Number of k-points (k x k x k)')
plt.ylabel('Total Energy (eV)')
plt.title('PBE KPOINTS test')
plt.savefig('vasp_kpts.png')
plt.close()
