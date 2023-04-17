# Path parser is applicated after the simulation, besides, it could be part of AutoSteper parameters.
# This is the pathway info for C60Cl6_1

import os

from autosteper.cage import Cage
from autosteper.parser import Path_Parser


a_cage = Cage(pristine_path=r'path/to/C60.xyz')
a_cage.set_workbase(root=r'path/to/C60_workbase')

path_para = {
    'step': 1,
    'start': 1,
    'q_add_num': 6,
    'q_path_rank': 5,
    'q_isomer_rank': 1,
    'log_low_e_num': 5,
}

a_path_parser = Path_Parser(path_para=path_para,
                            dump_root=r'path/to/pathway_dump_folder',
                            q_cage=a_cage)
a_path_parser.get_path_info()
