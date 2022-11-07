import os

from autosteper import AutoSteper

mach_para = {
        'batch_type': "x",
        'context_type': "SSHContext",
        'remote_root': 'xx/',
        'remote_profile': {
            "hostname": "xx",
            "username": "xx",
            "password": "xx",
            "port": 22,
            "timeout": 10
        }
    }

resrc_para = {
        'number_node': 6,
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
    'pristine_path': r'geom.xyz',
    'root': r'./',
    'gen_para': {'group': 'Cl',
                 'geom_mode': 'pre_defined',
                 'gen_core_path': r'xx\cagesearch.exe'
                 },
    'opt_mode': 'xtb',
    'opt_para': {
        'cmd_list': [r'xx/xtb', '--opt', 'tight', '--json'],
        'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
        'deal_wrong_mode': 'Report',
        'mach_para': mach_para,
        'resrc_para': resrc_para,
    },
    'run_para': {
        'start': 1,
        'stop': 4,
        'step': 1,
        'wht_list_para': {
            'mode': 'rank',
            'rank': 5
        }
    },
}

# raw folder is pre-defined
auto = AutoSteper(para=para)
auto.run()
