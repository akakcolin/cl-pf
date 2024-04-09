#!/usr/bin/env python


import os, sys, glob
import scipy
import numpy as np
import subprocess, fnmatch
import shutil
from pathlib import Path
from os import listdir
from os.path import isfile, join
from copy import deepcopy
import time
from optparse import OptionParser
from subprocess import Popen
import spglib
from ase.units import GPa, kJ
from ase.io import read, write
from ase import Atoms
from ase import units
from ase.calculators.dftb import Dftb

from ase.dft.kpoints import *
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry, is_subgroup
from ase.constraints import UnitCellFilter
from ase.optimize import BFGS
from ase.spacegroup import get_spacegroup
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo


usage = """%prog [options] xxx.cif"""
parser = OptionParser(usage)
(options, args) = parser.parse_args()
#atoms= read(args[0], format='vasp')
atoms= read(args[0], format='gen')


# We use an LJ calculator, and allow the cell and atomic positions to relax
b3_map = {'H': -0.1857, 'C': -0.1492, 'O': -0.1575, 'N': -0.1535,
        'Ca': -0.0340, 'K': -0.0339, 'Na': -0.0454, 'F': -0.1623,
        'Cl': -0.0697, 'Br': -0.0573, 'I': -0.0433, 'Mg': -0.02,
        'P': -0.14, 'S': -0.11, 'Zn': -0.03, 'B': -0.1008}
atoms_symbols = set([atom.symbol for atom in atoms])

b3_Ud = ["      %s = %8.4f" % (atom, b3_map[atom]) for atom in atoms_symbols]
dftb_d3="\n".join(b3_Ud)


#scaled_kpts=monkhorst_pack((2,2,2))
#scaled_kpts=get_monkhorst_pack_size_and_offset([[0, 0, 0]])
#abs_kpts=2*np.pi*np.dot(scaled_kpts, np.linalg.inv(atoms_prim.cell))

calc= Dftb(label='dftb',
            #kpts=mesh,
            Options_WriteResultsTag='Yes',
            Hamiltonian_SCC='Yes',
            Hamiltonian_SCCTolerance='1e-10',
            Hamiltonian_MaxSCCIterations =200,
            Hamiltonian_DampXH ='Yes',
            Hamiltonian_DampXHExponent =4.00,
            Hamiltonian_ThirdOrderFull='Yes',
            Hamiltonian_HubbardDerivs="{ \n" + dftb_d3+ "\n }",
            Hamiltonian_Dispersion ="SimpleDftD3 { \n a1=0.746 \n a2=4.191 \n s6=1.0 \n s8=3.209\n CutoffInter=94.868329 \n CoordinationNumber = exp {\n Cutoff=40 \n } \n }",
            #Driver_='ConjugateGradient',
            #Driver_MaxSteps=1000,
            #Driver_MaxForceComponent=0.001,
            #Driver_LatticeOpt='No',
            #Driver_FixCellOpt='Yes',
            #Driver_Pressure=1
            )
atoms.set_calculator(calc)
#print(atoms.get_potential_energy())
#atoms.set_constraint(FixSymmetry(atoms))
#atoms= UnitCellFilter(atoms)

#write('geo.gen', atoms, format='gen')
#dyn = BFGS(atoms, logfile='dyn.log')
potentialenergy = atoms.get_potential_energy()
vib = Vibrations(atoms)
vib.run()
res = vib.summary()
print(res)
vib_energies = vib.get_energies()

thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=atoms,
                        geometry='linear',
                        symmetrynumber=1, spin=0)
G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.)
harm= thermo.get_internal_energy(temperature=298.15)


print(G)
print(harm)
#dyn.run(fmax=0.001)
#print("Final Energy", atoms.get_potential_energy())

#write('geo_opted.gen', atoms, format='gen')
