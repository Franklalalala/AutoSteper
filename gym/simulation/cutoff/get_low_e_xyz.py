import numpy as np
import pandas as pd
from autosteper.tools import get_low_e_xyz
import os


para = {
    'mode': 'rank',
    'rank': 5
}

get_low_e_xyz(old_workbase=r'xx/AutoSteper/tests/test_step/geom',
              add_num=3,
              dump_folder=r'./dump',
              cutoff_para=para)

