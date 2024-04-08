#!/usr/bin/env python3
from ase.io import read, write
from ase.atoms import Atoms, Atom
import numpy as np
import spglib
from optparse import OptionParser


usage = """

# --------------------------------------------------------------------------------
# Convert conventional cell to primitive cell
# Aim: Read in CIF with ASE, use spglib to convert to an asymmetric unit cell
# Write asymmetric unit cell to xyz or CIF
# --------------------------------------------------------------------------------

%prog [option] XXX.cif """

parser = OptionParser(usage)
(options, args) = parser.parse_args()

# -------------------------------------------
# Functions
# -------------------------------------------

an_to_symbol = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
                10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K',
                20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
                30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y',
                40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In',
                50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr',
                60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm',
                70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
                80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac',
                90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',
                100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs',
                109: 'Mt',
                110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}

# https://github.com/atztogo/spglib/blob/master/python/examples/example.py
def show_symmetry(symmetry, n_symmetries=None):
    if n_symmetries == None:
        for i in range(symmetry['rotations'].shape[0]):
            print("  --------------- %4d ---------------" % (i + 1))
            rot = symmetry['rotations'][i]
            trans = symmetry['translations'][i]
            print("  rotation:")
            for x in rot:
                print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))

            print("  translation:")
            print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))
    else:
        for i in range(0, n_symmetries):
            print("  --------------- %4d ---------------" % (i + 1))
            rot = symmetry['rotations'][i]
            trans = symmetry['translations'][i]
            print("  rotation:")
            for x in rot:
                print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))

            print("  translation:")
            print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))


def show_lattice(lattice):
    print("Basis vectors:")
    for vec, axis in zip(lattice, ("a", "b", "c")):
        print("%s %10.5f %10.5f %10.5f" % (tuple(axis, ) + tuple(vec)))


def show_cell(lattice, positions, numbers):
    show_lattice(lattice)
    print("Atomic points:")
    for p, s in zip(positions, numbers):
        print("%2d %10.5f %10.5f %10.5f" % ((s,) + tuple(p)))


def asymmetric_cell_atom_indices(dataset):
    # atomic indices in supercell
    atom_indices = []
    for x in dataset['equivalent_atoms']:
        atom_indices.append(x)

    # Reduce to unique indices
    atom_indices = list(set(atom_indices))
    atom_indices.sort()
    return atom_indices


def spglib_to_ase(molecule, indices=None):
    basis = molecule[1]
    atomic_numbers = molecule[2]
    if indices == None:
        indices = range(0, len(atomic_numbers))

    #Have to store in Cartesian, not fractional (?)
    lattice = np.transpose(np.asarray(molecule[0]))
    print(lattice)

    ase_molecule = []
    for ia in indices:
        atomic_symbol = an_to_symbol[atomic_numbers[ia]]
        # Cartesian
        pos = np.matmul(lattice, np.asarray(basis[ia]))
        # Fractional
        #pos = basis[ia]
        print(ia, basis[ia])
        ase_molecule.append(Atom(atomic_symbol, pos))

    return Atoms(ase_molecule, cell=molecule[0], pbc=[1,1,1])

data_conventional = read(args[0], format="cif", store_tags=False)
print(vars(data_conventional))
print("")


data_primitive = read(args[0], format="cif", subtrans_included=False, primitive_cell=True)

# Primitive comes out looking erroneous, perhaps due to sublattice options being included in
# initial cif data
#write("ase-prim-aei.cif", data_primitive)
#write("ase-conv-aei.cif.cif", data_conventional)

# This hasn't made a super difference so continue with conventional

# Get ASE into spglib format. All tuples apart from atomic_numbers
lattice = []
for vector in data_conventional.cell:
    lattice.append(tuple(vector))

# ASE converts basis to cartesian, so convert back
inv_lattice = np.linalg.inv(np.transpose(np.asarray(data_conventional.cell)))

print(type(data_conventional))

basis = []
for pos in data_conventional.positions:
    frac_pos = np.matmul(inv_lattice, pos)
    basis.append(tuple(frac_pos))

atomic_numbers = []
for an in data_conventional.numbers:
    atomic_numbers.append(an)

print("Number of atoms", len(basis))

print("SGLIB data storage for a crystal:")
molecule = (lattice, basis, atomic_numbers)
# print(molecule)
print("  Spacegroup is %s." % spglib.get_spacegroup(molecule))

symmetry = spglib.get_symmetry(molecule)
# show_symmetry(symmetry)
print("  Number of symmetry operations is %d." % len(symmetry['rotations']))

# Only needs one rotation operation to determine the point group
# Not in my prefered notation. https://en.wikipedia.org/wiki/Crystallographic_point_group#Hermann–Mauguin_notation
print("  Pointgroup of aei is %s." %
      spglib.get_pointgroup(symmetry['rotations'])[0])

# From the international tables. Not sure how this function differs to get_symmetry
# BUT appear to need this apporach to get Wyckoff symbols/operations
dataset = spglib.get_symmetry_dataset(molecule)
print("  Spacegroup is %s (%d)." % (dataset['international'],
                                    dataset['number']))
print("  Pointgroup is %s." % (dataset['pointgroup']))
print("  Hall symbol is %s (%d)." % (dataset['hall'],
                                     dataset['hall_number']))
print("  Wyckoff letters are: ", dataset['wyckoffs'])

# print("  Mapping to equivalent atoms are: ")
# for i, x in enumerate(dataset['equivalent_atoms']):
#     print("  %d -> %d" % (i + 1, x + 1))




# Output from ASE in xyz
indices = asymmetric_cell_atom_indices(dataset)
ase_asymmetric_cell = spglib_to_ase(molecule, indices)
ase_asymmetric_cell.set_pbc((1, 1, 1))

#Need lattice vectors and add PBC
#Amm2 is base-centred orthorhombic
#https://en.wikipedia.org/wiki/Orthorhombic_crystal_system#Crystal_classes
# ----------------------------------------------------
# Reduce to primitive cell and find irreducible atoms
# ----------------------------------------------------
print(" Find primitive of conventional structure")
lattice, positions, numbers = spglib.find_primitive(molecule, symprec=1e-1)
show_cell(lattice, positions, numbers)
print("Number of atoms in primitive: ", len(numbers))
primitive_cell = (lattice, positions, numbers)

# https://atztogo.github.io/spglib/python-spglib.html#niggli-reduce
# The detailed control of standardization of unit cell can be done using standardize_cell.

ase_primitive_cell = spglib_to_ase(primitive_cell)
#write('aei_primtive_cell.xyz',ase_asymmetric_cell)


primitive_dataset = spglib.get_symmetry_dataset(primitive_cell)
print("  Spacegroup is %s (%d)." % (dataset['international'],
                                    dataset['number']))
print("  Pointgroup is %s." % (dataset['pointgroup']))
print("  Hall symbol is %s (%d)." % (dataset['hall'],
                                     dataset['hall_number']))
# print("  Wyckoff letters are: ", dataset['wyckoffs'])

print("  Mapping to equivalent atoms are: ")
for i, x in enumerate(primitive_dataset['equivalent_atoms']):
    print("  %d -> %d" % (i + 1, x + 1))

# Output from ASE in xyz
indices = asymmetric_cell_atom_indices(primitive_dataset)
ase_asymmetric_cell = spglib_to_ase(primitive_cell, indices)
write('asym_cell2.cif',ase_asymmetric_cell)
