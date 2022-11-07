import numpy as np
import pandas as pd
from autosteper.tools import get_low_e_ranks, get_low_e_xyz
import os

info = pd.read_pickle('passed_info.pickle')
e_arr = np.array(info['energy'])
#######################################################################################
para = {
    'mode': 'rank_or_value',
    'rank': 5,
    'value': 0.1
}
for a_rank in get_low_e_ranks(e_arr=e_arr, para=para, is_reverse=True):
    print(a_rank)
assert a_rank == len(e_arr) - 5
########################################################################################


#######################################################################################
para = {
    'mode': 'rank_and_value',
    'rank': 5,
    'value': 0.1
}
for a_rank in get_low_e_ranks(e_arr=e_arr, para=para, is_reverse=True):
    print(a_rank)
assert a_rank == len(e_arr) - 1
########################################################################################


#######################################################################################
para = {
    'mode': 'rank_and_value',
    'rank': 5,
    'value': 0.1
}
for a_rank in get_low_e_ranks(e_arr=e_arr, para=para):
    print(a_rank)
assert a_rank ==  0
########################################################################################


#######################################################################################
para = {
    'mode': 'rank_or_value',
    'rank': 5,
    'value': 0.1
}
for a_rank in get_low_e_ranks(e_arr=e_arr, para=para):
    print(a_rank)
assert a_rank == 4
########################################################################################

get_low_e_xyz(old_workbase=r'xx/AutoSteper/tests/test_step/geom',
              add_num=3,
              dump_folder=r'./dump',
              cutoff_para=para)

