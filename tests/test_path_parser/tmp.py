import os

import numpy as np
import pandas as pd


info = pd.read_pickle(r'all_info.pickle')
print(info)

# an_old = info.loc[2, 'e_list']
# print(an_old)
# an_old[2] = 1
# info.loc[2, 'e_list'] = an_old
# print(info.loc[2, 'e_list'])

e_arr = list(info['e_list'])
e_arr[0][0]=1
print(e_arr)
info['e_list'] = e_arr

