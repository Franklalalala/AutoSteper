from ase.units import Hartree, eV
from autosteper.parser import get_binding_e
from ase.io.gaussian import read_gaussian_out


H2 = read_gaussian_out(r'./H2.log')[-1]
H2_e = H2.get_total_energy()
Cage = read_gaussian_out(r'./cage.log')[-1]
Cage_e = Cage.get_total_energy()
get_binding_e(sorted_root=r'./sorted',
              addends_e=H2_e,
              cage_e=Cage_e)

