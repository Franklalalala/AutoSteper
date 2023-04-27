# Path parser is applicated after the simulation, besides, it could be part of AutoSteper parameters.
# This is the pathway info for C60Cl6_1

import os

from autosteper.parser import Path_Parser


refiner_para = {
    'has_parity': True,
    'group': 'Cl',
    'cmd_list': [r'path/to/xtb', '--opt', 'tight', '--json'],
    'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
    'deal_wrong_mode': 'Report',
    'mach_para': {
        'batch_type': "Torque",
        'context_type': "LocalContext",
        'remote_root': 'path/to/test_dpdispatcher',
        'remote_profile': None
    },
    'resrc_para': {
        'number_node': 1,
        'cpu_per_node': 6,
        'gpu_per_node': 0,
        'group_size': 1,
        'queue_name': "batch",
        'sub_batch_size': 10,
        'envs': {
            'dummy': 'dummy/path',
            # source your env here
            "OMP_STACKSIZE": "4G",
            "OMP_NUM_THREADS": "3,1",
            "OMP_MAX_ACTIVE_LEVELS": "1",
            "MKL_NUM_THREADS": "3"
        }
    },
}



a_path_parser = Path_Parser()
a_path_parser.update_path_info(old_workbase=r'path/to/sim_pathways',
                               new_workbase=r'path/to/updated_sim_pathways',
                               refine_top_k=4,
                               refine_mode='xtb',
                               refiner_para=refiner_para)
