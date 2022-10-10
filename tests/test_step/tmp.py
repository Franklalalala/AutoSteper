import numpy as np
import pandas as pd


info = pd.read_pickle('./passed_info.pickle')
aa = np.array(info['energy'])
print(aa)
ab = aa - min(aa)
for idx, a_e in enumerate(ab):
    if a_e > 1:
        break
print(idx)
print(ab)
range(len(ab))
rank_list = range