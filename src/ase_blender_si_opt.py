import bpy
import numpy as np
from ase.io import read, write
import blase as bl
from ase.calculators.dftb import Dftb
from ase.constraints import UnitCellFilter
from ase.optimize.sciopt import SciPyFminCG
from ase.optimize.bfgs import BFGS
from ase.dft.kpoints import *

from ase.calculators.loggingcalc import LoggingCalculator

atoms_prim = read("/test/si2_geoopt/POSCAR", format='vasp')


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


atoms_unit = atoms_prim.copy()
mesh, offset = genDftbKpoints(atoms_prim)
b3_map = {'H': -0.1857, 'C': -0.1492, 'O': -0.1575, 'N': -0.1535,
          'Ca': -0.0340, 'K': -0.0339, 'Na': -0.0454, 'F': -0.1623,
          'Cl': -0.0697, 'Br': -0.0573, 'I': -0.0433, 'Mg': -0.02,
          'P': -0.14, 'S': -0.11, 'Zn': -0.03, 'B': -0.1008}
atoms = set([atom.symbol for atom in atoms_prim])


b3_Ud = ["      %s = %8.4f" % (atom, b3_map[atom]) for atom in atoms]
dftb_d3 = "\n".join(b3_Ud)

# We use an LJ calculator, and allow the cell and atomic positions to relax
atoms_unit = atoms_prim.copy()

# scaled_kpts=get_monkhorst_pack_size_and_offset([[0, 0, 0]])
# abs_kpts=2*np.pi*np.dot(scaled_kpts, np.linalg.inv(atoms_prim.cell))
calc = Dftb(label='dftb',
            kpts=mesh,
            Hamiltonian_SCC='Yes',
            Hamiltonian_SCCTolerance='1e-5',
            Hamiltonian_MaxSCCIterations=200,
            Hamiltonian_DampXH='Yes',
            Hamiltonian_DampXHExponent=4.00,
            Hamiltonian_ThirdOrderFull='Yes',
            Hamiltonian_HubbardDerivs="{ \n" + dftb_d3 + "\n }",
            # Driver_='ConjugateGradient',
            # Driver_MaxSteps=2000,
            # Driver_MaxForceComponent=0.001,
            # Driver_LatticeOpt='Yes',
            # Driver_Isotropic='No',
            )

log_calc = LoggingCalculator(calc)
atoms_unit.set_calculator(log_calc)
ucf_unit = UnitCellFilter(atoms_unit)

atoms_unit.get_potential_energy()

dyn = BFGS(ucf_unit, trajectory='/test/si2_geoopt/dftb.traj',
           logfile='/test/si2_geoopt/dftb_unit.log')
# dyn = BFGS(ucf_unit, trajectory='dftb.traj', logfile='dftb_unit.log')
print("Unit Initial Energy", atoms_unit.get_potential_energy())
dyn.run(fmax=0.01)
print("Unit Final Energy", atoms_unit.get_potential_energy())
# write("dftb_unit.cif", atoms_unit, format='cif')

bl.clean_objects()


kwargs = {'show_unit_cell': 1,
          'radii': 0.2,
          'bond_cutoff': 1.0,
          'display': True
          }


kwargs2 = {'show_unit_cell': 1,
           'radii': 0.6,
           'bond_cutoff': 1.0,
           'display': True
           }

bl.view_blender(atoms_prim, **kwargs)
bl.view_blender(atoms_unit, **kwargs2)
