#!/usr/bin/env python

import os
import sys
import glob
import scipy
import numpy as np
import subprocess
import fnmatch
import shutil
from pathlib import Path
from os import listdir
from os.path import isfile, join
from copy import deepcopy
import time
from optparse import OptionParser
import spglib
from ase.units import GPa, kJ
from ase.io import read, write
from ase import Atoms
from ase import units
from ase.calculators.eam import EAM
from ase.calculators.aims import Aims
from ase.calculators.dftb import Dftb

from ase.dft.kpoints import *
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry, is_subgroup
from ase.constraints import UnitCellFilter
from ase.optimize import BFGS
from ase.spacegroup import get_spacegroup
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms


usage = """%prog [options] xxx.cif"""
parser = OptionParser(usage)
(options, args) = parser.parse_args()
# atoms= read(args[0], format='vasp')
atoms = read(args[0], format='gen')


def genDftbKpoints(lattice, h_max=15.):
    v = lattice.get_volume()
    va, vb, vc = lattice.get_cell()
    a, b, c, alpha, beta, gamma = lattice.get_cell_lengths_and_angles()
    ha = v / (np.sin(np.arccos(np.dot(vb, vc) / np.linalg.norm(vb) /
                               np.linalg.norm(vc)))
              ) / np.linalg.norm(vb) / np.linalg.norm(vc)
    hb = v / (np.sin(np.arccos(np.dot(va, vc) / np.linalg.norm(va) /
                               np.linalg.norm(vc)))
              ) / np.linalg.norm(va) / np.linalg.norm(vc)
    hc = v / (np.sin(np.arccos(np.dot(va, vb) / np.linalg.norm(va) /
                               np.linalg.norm(vb)))
              ) / np.linalg.norm(va) / np.linalg.norm(vb)
    # h_max = max([ha, hb, hc])
    x, y, z = map(lambda x: int(np.ceil(h_max/x)), [ha, hb, hc])
    nx = 0.0 if x % 2 == 1 else 0.5
    ny = 0.0 if y % 2 == 1 else 0.5
    nz = 0.0 if z % 2 == 1 else 0.5
    return (x, y, z), (nx, ny, nz)


# We use an LJ calculator, and allow the cell and atomic positions to relax
mesh, offset = genDftbKpoints(atoms)
b3_map = {'H': -0.1857, 'C': -0.1492, 'O': -0.1575, 'N': -0.1535,
          'Ca': -0.0340, 'K': -0.0339, 'Na': -0.0454, 'F': -0.1623,
          'Cl': -0.0697, 'Br': -0.0573, 'I': -0.0433, 'Mg': -0.02,
          'P': -0.14, 'S': -0.11, 'Zn': -0.03, 'B': -0.1008}
atoms_symbols = set([atom.symbol for atom in atoms])

b3_Ud = ["      %s = %8.4f" % (atom, b3_map[atom]) for atom in atoms_symbols]
dftb_d3 = "\n".join(b3_Ud)


# scaled_kpts=monkhorst_pack((2,2,2))
# scaled_kpts=get_monkhorst_pack_size_and_offset([[0, 0, 0]])
# abs_kpts=2*np.pi*np.dot(scaled_kpts, np.linalg.inv(atoms_prim.cell))

calc = Dftb(label='dftb',
            kpts=mesh,
            Options_WriteResultsTag='Yes',
            Hamiltonian_SCC='Yes',
            Hamiltonian_SCCTolerance='1e-10',
            Hamiltonian_MaxSCCIterations=200,
            Hamiltonian_DampXH='Yes',
            Hamiltonian_DampXHExponent=4.00,
            Hamiltonian_ThirdOrderFull='Yes',
            Hamiltonian_HubbardDerivs="{ \n" + dftb_d3 + "\n }",
            Hamiltonian_Dispersion="SimpleDftD3 { \n a1=0.746 \n a2=4.191 \n s6=1.0 \n s8=3.209\n CutoffInter=94.868329 \n CoordinationNumber = exp {\n Cutoff=40 \n } \n }",
            # Driver_='ConjugateGradient',
            # Driver_MaxSteps=1000,
            # Driver_MaxForceComponent=0.001,
            # Driver_LatticeOpt='No',
            # Driver_FixCellOpt='Yes',
            # Driver_Pressure=1
            )
atoms.set_calculator(calc)
print(atoms.get_potential_energy())


def build_dftb_qha(atoms, calc, step=0.02, shrink=-4, expand=5):
    run_dirs = []
    cif_names = []
    energies = []
    volumes = []
    for istep in range(shrink, expand):
        scale = np.power(1 + step * istep, 1.0 / 3.0)
        print(scale)
        shutil.rmtree('paralle_'+str(istep), ignore_errors=True)
        voldir = 'volume_'+str(istep)
        os.mkdir(voldir)
        atoms_it = atoms.copy()
        atoms_it.set_cell(atoms.get_cell() * scale, scale_atoms=True)
        atoms_it.set_calculator(calc)
        atoms_it.set_constraint(FixSymmetry(atoms_it))
        ucf_atoms = UnitCellFilter(atoms_it, constant_volume=True)
        dyn = BFGS(atoms_it)
        print("Initial Energy", atoms_it.get_potential_energy())
        print("Initial Volume", atoms_it.get_volume())
        dyn.run(fmax=0.0002)
        print("Final Energy", atoms_it.get_potential_energy())
        print("Final Volume", atoms_it.get_volume())
        print("initial symmetry at precision 1e-6")
        check_symmetry(atoms, 1.0e-6, verbose=True)
        print("sym symmetry after relaxation ")
        check_symmetry(atoms_it, 1.0e-6, verbose=True)
        write(voldir+'/dftb_opt_'+str(istep)+'.cif', atoms_it, format='cif')
        os.system('cp dftb_in.hsd_gamma ' + voldir)
        os.system('cp run_phonopy_displacment.py ' + voldir)
        forces = atoms_it.get_forces()
        energy = atoms_it.get_potential_energy()
        volume = atoms_it.get_volume()
        print(forces)
        print(energy)
        run_dirs.append(voldir)
        cif_names.append('dftb_opt_'+str(istep)+'.cif')
        energies.append(energy)
        volumes.append(volume)

    for i, dirs in enumerate(run_dirs):
        print(dirs)
        print(cif_names[i])
        os.chdir(dirs)
        print(os.getcwd())
        os.system('python run_phonopy_displacment.py   ' + cif_names[i])
        os.chdir('../')
    with open('e-v.dat', 'w') as f:
        f.write("# volume[A^3] - energy[eV]\n")
        for vol, en in sorted(zip(volumes, energies)):
            f.write("{:20.8f} {:20.8f}\n".format(vol, en))

    # Run volume-energy equation of state fitting
    os.system("phonopy-qha -b e-v.dat >fitting.stdout")
    with open("fitting.stdout") as f:
        fitting = yaml.safe_load(f.read())
    if not fitting:
        raise RuntimeError("Error occured while fitting e-v curve")
    if not (min(volumes) < fitting['Volume'] < max(volumes)):
        raise RuntimeError("Volume minima not in sampling range")


def serial_build_dftb_qha(atoms, calc, step=0.02, shrink=-4):
    run_dirs = []
    cif_names = []
    energies = []
    volumes = []
    scale = np.power(1 + step * shrink, 1.0 / 3.0)
    print(scale)
    shutil.rmtree('volume_'+str(shrink), ignore_errors=True)
    voldir = 'volume_'+str(shrink)
    os.mkdir(voldir)
    atoms_it = atoms.copy()
    atoms_it.set_cell(atoms.get_cell() * scale, scale_atoms=True)
    atoms_it.set_calculator(calc)
    atoms_it.set_constraint(FixSymmetry(atoms_it))
    ucf_atoms = UnitCellFilter(atoms_it, constant_volume=True)
    dyn = BFGS(atoms_it, logfile="dyn.log")
    dyn.run(fmax=0.0002)
    check_symmetry(atoms, 1.0e-6, verbose=True)
    print("sym symmetry after relaxation ")
    check_symmetry(atoms_it, 1.0e-6, verbose=True)
    write(voldir+'/dftb_opt_'+str(shrink)+'.cif', atoms_it, format='cif')
    os.system('cp dftb_in.hsd_gamma ' + voldir)
    os.system('cp run_phonopy_dftb.py ' + voldir)
    forces = atoms_it.get_forces()
    energy = atoms_it.get_potential_energy()
    volume = atoms_it.get_volume()
    print(forces)
    print(energy)
    run_dirs.append(voldir)
    cif_names.append('dftb_opt_'+str(shrink)+'.cif')
    energies.append(energy)
    volumes.append(volume)

    for i, dirs in enumerate(run_dirs):
        print(dirs)
        print(cif_names[i])
        os.chdir(dirs)
        print(os.getcwd())
        os.system('python run_phonopy_dftb.py   ' + cif_names[i])
        os.chdir('../')
    with open(str(shrink)+'_e-v.dat', 'w') as f:
        f.write("# volume[A^3] - energy[eV]\n")
        for vol, en in sorted(zip(volumes, energies)):
            f.write("{:20.8f} {:20.8f}\n".format(vol, en))

    # Run volume-energy equation of state fitting
    # os.system("phonopy-qha -b e-v.dat >fitting.stdout")
    # with open("fitting.stdout") as f:
    #    fitting = yaml.safe_load(f.read())
    # if not fitting:
    #    raise RuntimeError("Error occured while fitting e-v curve")
    # if not (min(volumes) < fitting['Volume'] < max(volumes)):
    #    raise RuntimeError("Volume minima not in sampling range")


def single_dftb_qha(atoms, calc):
    atoms_it = atoms.copy()
    atoms_it.set_calculator(calc)
    print(atoms_it.get_volume(), atoms_it.get_potential_energy())


def process_phonon_result(shrink=-4, expand=5):
    for istep in range(shrink, expand):
        voldir = 'volume_'+str(istep)
        os.chdir(voldir)
        os.system('phonopy -f disp-???/results.tag --dftb+')
        os.system('phonopy -t mesh.conf --dftb+')
        os.system(
            'cp thermal_properties.yaml ../thermal_properties.yaml_'+str(istep))
        os.chdir('../')


# build_dftb_qha(atoms, calc, step=0.01, shrink=-4, expand=5)
serial_build_dftb_qha(atoms, calc, step=0.01, shrink=1)

# single_dftb_qha(atoms, calc)
# process_phonon_result(shrink=-4, expand=6)
