import os

from autosteper import AutoSteper


para = {
    'pristine_path': r'geom.xyz',
    'root': r'./',
    'gen_para': {'group': 'CF3',
                 'geom_mode': 'pre_defined',
                 'gen_core_path': r'xx\cagesearch.exe'
                 },
    'opt_mode': 'xtb',
    'opt_para': {
        'cmd_list': [r'xx/xtb', '--opt', 'tight', '--json', '--cycles 10'],
        'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
        'deal_wrong_mode': 'Report',
        'mach_para': {
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
        },
        'resrc_para': {
            'number_node': 6,
            'cpu_per_node': 6,
            'gpu_per_node': 0,
            'group_size': 10,
            'queue_name': "batch",
            'sub_batch_size': 50
        },
    },
    'random_para': {
        'addon_list': [30, 12],
        'random_num': 8,
        'try_times': 3
    }
}



# raw folder is pre-defined
auto = AutoSteper(para=para)
auto.random()
