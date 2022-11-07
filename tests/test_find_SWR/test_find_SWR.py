import os

from autosteper.parser import find_SWR


find_SWR(q_sorted_root=r'./query_workbase',
         tgt_sorted_root=r'./target_workbase',
         swr_dump_path=r'./SWR_dump',
         step=2,
         is_unique=False,
         is_low_e=False)




