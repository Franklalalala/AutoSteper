import os

import numpy as np
import pandas as pd
from autosteper.tools import get_low_e_ranks

info = pd.read_pickle('passed_info.pickle')
e_arr = np.array(info['energy'])  # or e_arr = info['energy']

#######################################################################################
para = {
    'mode': 'value',
    'value': 0.1
}
for a_rank in get_low_e_ranks(e_arr=e_arr, para=para):
    print(a_rank)
assert a_rank == 0

#######################################################################################
para = {
    'mode': 'rank',
    'rank': 5
}
for a_rank in get_low_e_ranks(e_arr=e_arr, para=para):
    print(a_rank)
assert a_rank + 1 == 5

#################################### lower boundary of rank, value ########################################
# met any one of the cutoff will be good
para = {
    'mode': 'rank_and_value',
    'rank': 5,
    'value': 0.1
}
for a_rank in get_low_e_ranks(e_arr=e_arr, para=para):
    print(a_rank)
assert a_rank == 0

################################ upper boundary of rank, value ###############################################
# both cutoff needs to be met
para = {
    'mode': 'rank_or_value',
    'rank': 5,
    'value': 0.1
}
for a_rank in get_low_e_ranks(e_arr=e_arr, para=para):
    print(a_rank)
assert a_rank + 1 == 5

################################ reverse cutoff for black list ###############################################
para = {
    'mode': 'rank',
    'rank': 2,
}
for a_rank in get_low_e_ranks(e_arr=e_arr, para=para, is_reverse=True):
    print(a_rank)
assert len(e_arr) - para['rank'] == a_rank
