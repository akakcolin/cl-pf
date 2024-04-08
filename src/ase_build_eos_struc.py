
from ase import *
from ase.io import read,write
from ase.units import kJ, _e
import numpy as np
from math import sqrt
# setting vasp parameter
refcell = np.array([[1.0, 0.0, 0.0],
		 [0.0, 1.0, 0.0],
		[0.0, 0.0, 1.0]])

# read initial POSCAR

f = file('file')
line = f.readline().strip('\n')
filename=line
print filename
try:
	atoms = read(filename)
	print atoms
	init_cell=atoms.get_cell()
	print init_cell
except:
	print "IOError: Could not determine chemical symbols"
	pass
label=0
for x in np.linspace(0.85, 1.20, 10):
	atoms.set_cell(init_cell*x, scale_atoms=True)
	write(str(label)+'_'+filename, atoms,format="vasp", direct=True, sort=True)
	label=label+1
