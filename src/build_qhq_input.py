import numpy as np
import ase
from ase.visualize import view
from ase.spacegroup import crystal, Spacegroup
from ase.build import molecule
from ase.io import read, write

def is_symmorphic(trans):
    for t in trans:
        if(sum(t)):
            return 0
    return 1

h2o=molecule('H2O')
h2o_positions=h2o.get_scaled_positions()
h2o_symbols=h2o.get_chemical_symbols()
mol_unit=[]
shift_vector=[0.28, 0.3, 0.32]

new_positions=[ pos+shift_vector for pos in h2o_positions]

h2o_asym=zip(h2o_symbols, new_positions)
for element in h2o_asym:
    t_atoms = ase.Atom(element[0], element[1])
    mol_unit.append(t_atoms)

print(mol_unit)

a = b = c=4
alpha = beta = 90
gamma = 90
cells=[]

def print_head(is_symmorphic):
    if(is_symmorphic):
        return "1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1"
    else:
        return "1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 0"

def print_primitive_cell(primiteve_cell):
    string_cell=""
    for i in primiteve_cell:
       # print("{0} {1} {2}".format(i[0], i[1], i[2]))
        string_cell += (str(i[0]) + " ")
        string_cell += (str(i[1]) + " ")
        string_cell += (str(i[2]) + " ")
        string_cell += "\n"
    return string_cell[0:-1]

def get_each_number(chemical_symbols):
    b={}
    for i in chemical_symbols:
        if i not in b:
            b[i]=1
        else:
            b[i]+=1
    return b

def print_translation(trans):
    for i in trans:
        print("0")
        print("{0} {1} {2}".format(i[0], i[1], i[2]))

def get_point_name(op):

    point_numbers = [0,0,
            2,
            3,3,3,
            4,4,4,4,
            5,5,5,5,5,
            6,6,6,6,6,6,6,6,6,
            7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
            8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
            9,9,9,9,9,9,
            10]

    return point_numbers[op]

for i in range(30):
    filename="h2o_" + str(i+2)
    spacegroup = Spacegroup(i+2)
    f=open(filename+".dat", "w")
    dftb_f=open(filename+".gen", "w")

    rot, trans=spacegroup.get_op()
    #print("translations:")
    #print(trans)
    #print("is symmorphic :",is_symmorphic(trans))
    print(print_head(is_symmorphic(trans)), file=f)
    #print("scaled_primitive_cell")
    #print(spacegroup.scaled_primitive_cell)
    print(print_primitive_cell(spacegroup.scaled_primitive_cell), file=f)
    print(get_point_name(spacegroup.no), file=f)
    #print('Space roup', spacegroup.no, spacegroup.symbol, file=f)
    print("2\n1\n0", file=f)  #HO
    #print("3\n1\n1\n0", file=f)  #HO
    #print("1\n1", file=f)   # N

    cell = crystal(mol_unit, spacegroup=spacegroup, cellpar=[a, b, c, alpha, beta, gamma], primitive_cell=True)
    #print(len(cell.get_scaled_positions()))
    #print(cell.get_chemical_symbols(), file=f)
    num_b=get_each_number(cell.get_chemical_symbols())
    #print(num_b['C'], file=f)
    print(num_b['O'], file=f)
    print(num_b['H'], file=f)
    cell_positions=cell.get_scaled_positions()
    a=1
    chemical_symbols=cell.get_chemical_symbols()
    print(" {0} F".format(len(cell.get_chemical_symbols())), file=dftb_f)
    print(" O H ", file=dftb_f)
    for i in range(len(cell.get_scaled_positions())):
        print("0", file=f)
        print("{0} {1} {2}".format(cell_positions[i][0], cell_positions[i][1], cell_positions[i][2]), file=f)
        if chemical_symbols[i] == 'O' :
            b=1
        elif chemical_symbols[i] == 'H':
            b=2
        else:
            b=3
        print("{0} {1} {2} {3} {4}".format(a, b, cell_positions[i][0], cell_positions[i][1], cell_positions[i][2]), file=dftb_f)
        a=a+1
    print("0.0000000 0.0000000 0.00000", file=dftb_f)
    print(print_primitive_cell(spacegroup.scaled_primitive_cell), file=dftb_f)
    dftb_f.close()

    if(np.count_nonzero(trans)):
          for i in trans:
              print("0", file=f)
              print("{0} {1} {2}".format(i[0], i[1], i[2]), file=f)


    print("1\n0\n0", file=f)
    print("0.0 0.0 0.0", file=f)
    print("0\n1\n0\n0.0 0.0 0.0 \n0", file=f)
    f.close()

    #print(cell.get_scaled_positions())

    #write(filename+".cif", cell, format='cif')
    #write(filename+".vasp", cell, format='vasp', direct=True, sort=False)

    #cells.append(cell)
#view(cell)
