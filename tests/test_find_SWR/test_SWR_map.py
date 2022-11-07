import os

from ase.io.gaussian import read_gaussian_out
from autosteper.tools import map_SWR


q_atoms = read_gaussian_out(fd=r'./query_workbase/log/0_addons_1.log')[-1]
tgt_atoms = read_gaussian_out(fd=r'./target_workbase/log/0_addons_1.log')[-1]
a_swr = map_SWR(q_atoms=q_atoms, tgt_atoms=tgt_atoms)
print(a_swr)