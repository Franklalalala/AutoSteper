import os

from autosteper.parser import find_SWR, count_SWR

# set is_low_e on to enable an energy criterion
# set is_unique on to keep only one swr product
find_SWR(q_sorted_root=r'./11_C84',
         tgt_sorted_root=r'./12_C84',
         swr_dump_path=r'./11_to_12',
         step=2,
         is_unique=True,
         is_low_e=True)

find_SWR(q_sorted_root=r'./12_C84',
         tgt_sorted_root=r'./11_C84',
         swr_dump_path=r'./12_to_11',
         step=2,
         is_unique=True,
         is_low_e=True)
