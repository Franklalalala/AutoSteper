import numpy as np
import pandas as pd
from ase.units import Hartree


info = pd.read_pickle('./deep_yes_info.pickle')
aa = np.array(info['energy'])*Hartree
print(aa[:5])
# ab = aa - min(aa)
# for idx, a_e in enumerate(ab):
#     if a_e > 1:
#         break
# print(idx)
# print(ab)
# range(len(ab))
# rank_list = range
