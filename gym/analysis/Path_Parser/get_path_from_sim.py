# Path parser is applicated after the simulation, besides, it could be part of AutoSteper parameters.
# This is the pathway info for C60Cl6_1

import os

from autosteper.parser import Path_Parser


a_path_parser = Path_Parser()
a_path_parser.get_path_from_sim(dump_root=r'path/to/sim_pathways',
                                pristine_cage_path=r'path/to/C60.xyz',
                                sim_workbase=r'path/to/C60',
                                sim_step=1, sim_start=1, q_add_num=6, q_path_rank=5, q_isomer_rank=3,
                                is_mix=True)
