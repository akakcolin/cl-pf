#!/usr/bin/env python
import os
import sys

import numpy as np
from optparse import OptionParser

from ase.io import read, write
from ase import Atoms

from ase.calculators.dftb import Dftb
from ase.calculators.vasp import Vasp
from ase.spacegroup.symmetrize import FixSymmetry
from ase.constraints import UnitCellFilter
from ase.optimize import BFGS, FIRE, MDMin
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG

usage = """%prog [options] xxx.cif"""
parser = OptionParser(usage)
(options, args) = parser.parse_args()
#atoms= read(args[0], format='vasp')
atoms= read(args[0], format='cif')

kpar = 1
ncore=24

calc = Vasp(setups='recommended',
            xc='PBE',           #PAW potentials [PBE or LDA]
            prec='normal',      #specifies the "precision"-mode [Accurate strictly avoids any aliasing or wrap around errors]
            ediff=1e-6,         #Convergence criterion [eV]
            nsw=0,              #maximum number of ionic steps
            ibrion=-1,          # [-1: no update, 0: molecular dynamics, 1: ionic relaxation (RMM-DIIS) ...]
            nbands=30,          #number of orbitals [default: NELECT/2+NIONS/2	non-spinpolarized]
            lcharg='.TRUE.',    #determines whether the charge densities (files CHGCAR and CHG) are written
            icharg = 1,         #determines how VASP constructs the initial charge density
            kpts=[5, 5, 5],     #specifies the Monkhorst-Pack grid
            ismear=1,           #Methfessel Paxton energy smearing order

            sigma=0.2)          #Smearing energy [eV]

calc.set(atoms=atoms)           #Set ASE Atoms object
calc.set(kpar=int(kpar), ncore=int(ncore)) #Set kpar and ncore
atoms.set_calculator(calc)
print(atoms.get_potential_energy())
