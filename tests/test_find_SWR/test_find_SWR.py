import os

from autosteper.parser import find_SWR


find_SWR(q_sorted_root=r'F:\final_paper_data\SWR\workbase\16',
         tgt_sorted_root=r'F:\final_paper_data\SWR\workbase\11',
         swr_dump_path=r'F:\AutoSteper\tests\test_find_SWR\new_SWR_dump',
         step=2,
         is_unique=False,
         is_low_e=False)




