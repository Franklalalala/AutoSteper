import os

from autosteper.parser import find_SWR


find_SWR(q_logs=r'F:\AutoSteper\tests\test_find_SWR\query_workbase\disordered_logs',
         tgt_logs=r'F:\AutoSteper\tests\test_find_SWR\target_workbase\disordered_logs',
         log_mode='gauss',
         q_workbase=r'./query_workbase',
         tgt_workbase=r'./target_workbase',
         swr_dump_path=r'./SWR_dump',
         step=2,
         is_low_e=True,
         is_unique=True,
         file_mid_name='dft')

