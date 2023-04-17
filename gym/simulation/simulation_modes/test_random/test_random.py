# Note: here performed random test on two C60(CF3)m, that is m = 12 and m = 24.
# For m = 12, all 8 randomly generated isomers have passed topological check,
# as there are 8 xyz files in ./C60/12addons/cooked.
# For m = 24, only 4 isomers passed topological check, see ./C60/24addons/cooked
# However, there are only 2 isomers in the ./C60/24addons/failed_job_paths
# The rest two isomers went "wrong". No output files received from remote.
# See the ./C60/24addons/status_info.pickle

import os

from autosteper import AutoSteper

mach_para = {
    'batch_type': "Torque",
    'context_type': "LocalContext",
    'remote_root': 'path/to/test_dpdispatcher',
    'remote_profile': None
}

resrc_para = {
    'number_node': 'node05',
    'cpu_per_node': 6,
    'gpu_per_node': 0,
    'group_size': 1,
    'queue_name': "batch",
    'envs': {
        "OMP_STACKSIZE": "4G",
        "OMP_NUM_THREADS": "3,1",
        "OMP_MAX_ACTIVE_LEVELS": "1",
        "MKL_NUM_THREADS": "3"
    },
    'sub_batch_size': 50
}


para = {
    'pristine_path': r'C60.xyz',
    'root': r'./',
    'gen_para': {'group': 'CF3',
                 'geom_mode': 'pre_defined',
                 'gen_core_path': r'path/to/cagesearch'
                 },
    'opt_mode': 'xtb',
    'opt_para': {
        'cmd_list': [r'path/to/xtb', '--opt', 'tight', '--json'],
        'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
        'deal_wrong_mode': 'Ignore',
        'mach_para': mach_para,
        'resrc_para': resrc_para,
    },
    'random_para': {
        'addon_list': [12, 24],
        'random_num': 8,
        'try_times': 3
    }
}

auto = AutoSteper(para=para)
auto.random()
