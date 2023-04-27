import os

from autosteper import AutoSteper

mach_para = {
    'batch_type': "Torque",
    'context_type': "LocalContext",
    'remote_root': 'x/to/test_dpdispatcher',
    'remote_profile': None
}

resrc_para = {
    'number_node': 1,
    'cpu_per_node': 6,
    'gpu_per_node': 0,
    'group_size': 10,
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
    'gen_para': {'group': 'Cl',
                 'geom_mode': 'pre_defined',
                 'gen_core_path': r'path/to/cagesearch'
                 },
    'opt_mode': 'xtb',
    'opt_para': {
        'has_parity': True,
        'cmd_list': [r'path/to/xtb', '--opt', 'tight', '--json'],
        'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
        'deal_wrong_mode': 'Report',
        'mach_para': mach_para,
        'resrc_para': resrc_para,
    },
    'run_para': {
        'start': 1,
        'stop': 3,
        'step': 1,
        'wht_list_para': {
            'mode': 'rank',
            'rank': 3
        }
    },
}

auto = AutoSteper(para=para)
auto.run()
