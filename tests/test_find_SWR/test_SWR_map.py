import os

from ase.io.gaussian import read_gaussian_out
from autosteper.tools import map_SWR


q_atoms = read_gaussian_out(fd=r'F:\AutoSteper\tests\test_find_SWR\query_workbase\disordered_logs\0_addons_dft_0.log')[-1]
tgt_atoms = read_gaussian_out(fd=r'F:\AutoSteper\tests\test_find_SWR\target_workbase\disordered_logs\0_addons_dft_0.log')[-1]
a_swr = map_SWR(q_atoms=q_atoms, tgt_atoms=tgt_atoms)
print(a_swr)