# Path parser is applicated after the simulation, besides, it could be part of AutoSteper parameters.
# This is the pathway info for C60Cl6_1

import os

from autosteper.parser import Path_Parser


a_path_parser = Path_Parser()
a_path_parser.re_generate_pathways(old_workbase=r'path/to/updated_sim_pathways',
                                   new_workbase=r'path/to/regenerated_pathways',
                                   step=1, last_log_mode='xtb', keep_top_k_pathway=3, group='Cl')
