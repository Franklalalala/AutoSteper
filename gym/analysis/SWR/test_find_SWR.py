import os

from autosteper.parser import find_SWR


# set is_low_e on to enable an energy criterion
# set is_unique on to keep only one swr product
find_SWR(q_sorted_root=r'./cook_disordered_results/11/sorted',
         tgt_sorted_root=r'./cook_disordered_results/12/sorted',
         swr_dump_path=r'./q_11_to_tgt_12',
         is_unique=True,
         is_low_e=True)

