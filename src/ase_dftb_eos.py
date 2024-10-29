from ase import Atoms

from ase.io.trajectory import Trajectory

from ase.eos import EquationOfState
from ase.io import read
from ase.units import kJ

from ase.io.xyz import read_xyz,  write_xyz
from ase.io import read, write

import time
import numpy as np
import os, shutil
from ase.io import Trajectory, read, write, iread

from ase.atoms import *
from ase.calculators.dftb import Dftb
from ase.calculators.singlepoint import SinglePointCalculator

from ase.io import extxyz
import sys
import os
import numpy as np
from ase.units import Bohr, Hartree

from tqdm import tqdm
import numpy as np
import argparse
import itertools
import sys
import os
import os.path
import tempfile
import shutil
import ase.io
import ase.atoms

os.environ['DFTB_COMMAND'] = '/usr/local/bin/dftb+'
base_dir = os.getcwd()
os.chdir(base_dir)

print('You are here:', base_dir)
# set up model for the DFTB+ computations using ASE


def genDftbKpoints(lattice, h_max=15.):
    v = lattice.get_volume()
    va, vb, vc = lattice.get_cell()
    a, b, c, alpha, beta, gamma = lattice.get_cell_lengths_and_angles()
    ha = v / (np.sin(np.arccos(np.dot(vb, vc) / np.linalg.norm(vb) /
              np.linalg.norm(vc)))) / np.linalg.norm(vb) / np.linalg.norm(vc)
    hb = v / (np.sin(np.arccos(np.dot(va, vc) / np.linalg.norm(va) /
              np.linalg.norm(vc)))) / np.linalg.norm(va) / np.linalg.norm(vc)
    hc = v / (np.sin(np.arccos(np.dot(va, vb) / np.linalg.norm(va) /
              np.linalg.norm(vb)))) / np.linalg.norm(va) / np.linalg.norm(vb)
    # h_max = max([ha, hb, hc])
    x, y, z = map(lambda x: int(np.ceil(h_max/x)), [ha, hb, hc])
    nx = 0.0 if x % 2 == 1 else 0.5
    ny = 0.0 if y % 2 == 1 else 0.5
    nz = 0.0 if z % 2 == 1 else 0.5
    return (x, y, z), (nx, ny, nz)


def Kgrid(struc, Kresol=0.10, dimension=3):
    """
    Assign kpoints based on the lattice
    """
    a, b, c, alpha, beta, gamma = struc.get_cell_lengths_and_angles()
    vol = struc.get_volume()
    dist = np.zeros(3)
    dist[2] = np.abs(vol/(a*b*np.sin(np.radians(gamma))))
    dist[1] = np.abs(vol/(a*c*np.sin(np.radians(beta))))
    dist[0] = np.abs(vol/(b*c*np.sin(np.radians(alpha))))
    Kpoints = np.ceil(1./(dist*Kresol))
    if dimension == 2:
        Kpoints[-1] = 1
    return Kpoints.astype(int)

Elements_Angulars = {
    "Br": "d", "Mg": "p",
    "C":  "p", "N":  "p",
    "Ca": "p", "Na": "p",
    "Cl": "d", "O":  "p",
    "F":  "p", "P":  "p",
    "H":  "s", "S":  "d",
    "I":  "d", "Zn": "d",
    "K":  "p", "B": "p",
    "Cu": "d", "Pd": "d",
    "Si": "p", "Li": "s",
    "Hf": "f"
}


def gen_max(angulars):
    hamil_max = ''
    for i, b in angulars:
        hamil_max = hamil_max + "             "+i + "=\"" + b + "\"\n"
    return hamil_max

def setup_dftb_calc(struct, charge=0, **kwargs):
    #print('setting up DFTBplus calculator')
    # link to your DFTB binary. Modify if needed.

    os.environ["DFTB_PREFIX"] = "./fit_skf/"
    mesh, offset = genDftbKpoints(struct)

    atoms_symbols = set([atom.symbol for atom in struct])

    angulars = [(at, Elements_Angulars[at]) for at in atoms_symbols]
    ang_str = gen_max(angulars)

    DFTB_calc = Dftb(label='dftb', kpts=mesh,
                     Options_WriteResultsTag='Yes',
                     Hamiltonian_MaxSCCIterations='500',
                     Hamiltonian_SCC='Yes',
                     Hamiltonian_SCCTolerance='1e-6',

                     Hamiltonian_MaxAngularMomentum="{ \n" + \
                     ang_str + "     \n}",
                     Hamiltonian_Filling= "MethfesselPaxton { Temperature=0.01 }",
                     )

    if charge > 0:
        DFTB_calc.set(Hamiltonian_Charge=str(charge))

    calc = DFTB_calc
    return calc

usage="""test """
parser = argparse.ArgumentParser(prog="build_dataset.py", usage=usage)
parser.add_argument("inp", type=str, default="geometries.xyz",
 help="XYZ file containing geometries [default: geometries.xyz]")

args = parser.parse_args()
atoms = ase.io.read(args.inp)

calculator = setup_dftb_calc(atoms)
atoms.calc = calculator
cell = atoms.get_cell()

traj = Trajectory('dftb_eos.traj', 'w')
volumes =[]
energies =[]
for x in np.linspace(0.80, 1.20, 10):
    atoms.set_cell(cell * x, scale_atoms=True)
    eng = atoms.get_potential_energy()
    energies.append(eng)
    vol = atoms.get_volume()
    volumes.append(vol)
    traj.write(atoms)
print(volumes)
print(energies)
# Extract volumes and energies:

eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
eos.plot('dftb_eos.png')
