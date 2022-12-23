import os

from autosteper import AutoSteper
from autosteper.cage import Cage
from autosteper.parser import Path_Parser


a_cage = Cage(pristine_path=r'./geom.xyz')
a_cage.set_workbase(root=r'xx')

path_para = {
    'step': 1,
    'start': 1,
    'q_add_num': 4,
    'q_path_rank': 10,
    'q_isomer_rank': 1,
    'log_low_e_num': 10,

}

a_path_parser = Path_Parser(path_para=path_para,
                            dump_root=r'xx',
                            q_cage=a_cage)
a_path_parser.get_path_info()
