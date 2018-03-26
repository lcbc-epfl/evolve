from ase import Atoms
from ase.calculators.amber import Amber

atoms = Atoms('OH2OH2',
              [[-0.956, -0.121, 0],
               [-1.308, 0.770, 0],
               [0.000, 0.000, 0],
               [3.903, 0.000, 0],
               [4.215, -0.497, -0.759],
               [4.215, -0.497, 0.759]])

calc = Amber(amber_exe='sander -O',
             infile='mm.in',
             outfile='mm.out',
             topologyfile='2H2O.top',
             incoordfile='mm.crd')
calc.write_coordinates(atoms, 'mm.crd')
atoms.set_calculator(calc)
f = atoms.get_forces()

print f

print atoms.get_kinetic_energy()

print atoms.get_potential_energy()