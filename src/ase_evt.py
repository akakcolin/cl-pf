"""Demonstrates molecular dynamics with constant energy."""

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton,fire,gpmin,mdmin,LBFGS
from ase.io import write,read
from ase.build import molecule
from rdkit import Chem
from rdkit.Chem import AllChem
from ase.optimize.minimahopping import MinimaHopping
from ase import *
from ase.optimize.basin import BasinHopping
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase import units

# Set the momenta corresponding to T=300K
atoms = read('trail.mol')
atoms.set_calculator(Dftb(label='h2o',
                         atoms=atoms,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_C='"p"',
                         Hamiltonian_MaxAngularMomentum_N='"p"',
                         Hamiltonian_MaxAngularMomentum_O='"p"',
                         Hamiltonian_MaxAngularMomentum_H='"s"',
                         ))
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
# We want to run MD with constant energy using the Langevin algorithm
# with a time step of 5 fs, the temperature T and the friction
# coefficient to 0.02 atomic units.
T = 300
dyn = Langevin(atoms, 1 * units.fs, T * units.kB, 0.002)


def printenergy(a=atoms):  # store a reference to atoms in the definition.
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

dyn.attach(printenergy, interval=50)

# We also want to save the positions of all atoms after every 100th time step.
traj = Trajectory('moldyn3.traj', 'w', atoms)
dyn.attach(traj.write, interval=50)

# Now run the dynamics
printenergy()
dyn.run(5000)
