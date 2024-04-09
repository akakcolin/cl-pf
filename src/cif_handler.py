#!/usr/bin/python3
import os
from ase.io import read, write
import numpy as np

import itertools
import re
import numpy as np
import scipy.stats as stats
from collections import OrderedDict
import argparse
import pathlib
import tempfile
import os
import subprocess
import copy
from pymatgen.core import Structure
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from copy import deepcopy

general_format='-i {cif} -s {cell} -v 1 -q -m -n l{smp} -n h{smp} -n r{smp}'

def run_supercell(cmd, params, outdir):
    fullcmd="{cmd} {params} -o {out}".format(cmd=cmd, params=params, out=outdir / "electro")
    print(fullcmd)
    result = subprocess.run(fullcmd, capture_output=True, text=True, shell=True)

    if result.returncode != 0:
        print(result.stderr)
        return None, None

    tot_syms = re.search("([0-9]+) +symmetry operation found for supercell.*", result.stdout)
    tot_comb = re.search("Combinations after merge: ([0-9]+)", result.stdout)
    return outdir, int(tot_syms.group(1)), int(tot_comb.group(1))

def get_cif_files(cifpath, skip_files=None):
    """
        Get the list of CIF files
        Args:
                mofpath (string): directory to CIF files
                skip_files (list): list of file to ignore
        Returns:
                sorted_cifs (list): alphabetized list of CIF files
        """
    cif_files = []
    if skip_files is None:
        skip_files = []
    for filename in os.listdir(cifpath):
        filename = filename.strip()
        if '.cif' in filename or 'POSCAR_' in filename:
            if '.cif' in filename:
                refcode = filename.split('.cif')[0]
            elif 'POSCAR_' in filename:
                refcode = filename.split('POSCAR_')[1]
            if refcode not in skip_files:
                cif_files.append(cifpath+'/'+filename)
            else:
                print('Skipping '+refcode)


    sorted_cifs = sorted(cif_files)
    return sorted_cifs

def cif_to_aseobj(filepath):
    """
        Convert file to ASE Atoms object
        Args:
                filepath (string): full path to structure file
        Returns:
                sorted_cifs (list): alphabetized list of CIF files
        """
    tol = 0.8
    aseobj = read(filepath)
    d = aseobj.get_all_distances()
    min_val = np.min(d[d>0])

    if min_val < tol:
        print('WARNING: Atoms overlap by '+str(min_val))
    return aseobj

def substitute_atom_by_atom_symbol(atoms, atom_idxs, new_atom_symbol):
    new_atoms = atoms.copy()
    print(atoms.symbols)
    for idx in atom_idxs:
        new_atoms[idx].symbol = new_atom_symbol
        print("old atom's symbol")
        print(atoms[idx].symbol)
        print(new_atoms[idx].symbol)
    return new_atoms

def get_substitue_index_by_symbol(atoms, atom_symbol):

    atoms_symbols = atoms.get_chemical_symbols()
    indices=[]
    for i in range(len(atoms_symbols)):
        if atom_symbol == atoms_symbols[i]:
            indices.append(i)
    return indices

def substitute_atom_by_molecular(atoms, atom_idxs, molecule):
    new_atoms = atoms.copy()
    cell = atoms.cell

    positions=atoms.get_scaled_positions()
    print(atoms.symbols)
    cents=[]
    for idx in atom_idxs:
        cents.append(positions[idx])
    print(cents)

    for ref in cents:
        new_mol = deepcopy(molecule)
        for atom in new_mol:
            atom.position = atom.position + ref
            print(atom.position)
        new_atoms.extend(new_mol)
        del new_atoms
    atoms_symbols = new_atoms.get_chemical_symbols()

    for i in range(len(atoms_symbols)):
        if atom_symbol != atoms_symbols[i]:
            final_atoms.append(new_atoms[i])

    #for idx in atom_idxs:
    #    del new_atoms[idx]
    return new_atoms

def substitute_atom_by_molecular2(atoms, atom_idx):
    new_atoms = atoms
    new_atoms.set_atom_number
    new_atoms.replace("I", "Br", num_sub)
    new_atoms.to_vasp_poscar("POSCAR_substituted"+str(i)+".vasp")

    coords = molecular.center_coord


if __name__ == '__main__':
#    parser = argparse.ArgumentParser(description='Create tests.')
#    parser.add_argument('--supercell-cmd', type=str, default='supercell', help='Path to checked supercell program')
#    parser.add_argument('--data-folder', type=pathlib.Path, default=pathlib.Path("./"), help='Path with data files.')
#    parser.add_argument('--tmp-folder', type=pathlib.Path, default=pathlib.Path(tempfile.gettempdir()), help='Path to store temporary files')
#    parser.add_argument('-v', '--verbose', action='store_true')

#    args = parser.parse_args()
#    print(args)
#    if args:
#        main(args)
    ciffiles = get_cif_files('test')
    for i in ciffiles:
        print(i)

    first_struc=read(ciffiles[0])

    from ase.build import molecule
    ch4 = molecule('CH4')


    indice_ca=get_substitue_index_by_symbol(first_struc, 'Ca')
    print(indice_ca)
    #new_struct = substitute_atom_by_atom_symbol(first_struc, 1, 'Cs')
    new_struct=substitute_atom_by_molecular(first_struc, indice_ca, ch4)

    write('new_struct.cif', new_struct)
