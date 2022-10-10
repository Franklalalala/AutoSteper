import os

import numpy as np
import pandas as pd


info = pd.read_pickle(r'passed_info.pickle')
e_list = info['energy']
e_arr = np.array(e_list)
e_arr = e_arr - min(e_arr)
for idx, value in enumerate(e_arr):
    if value > 1:
        break
print(idx)
print(e_arr[idx])
print(e_arr[idx-1])