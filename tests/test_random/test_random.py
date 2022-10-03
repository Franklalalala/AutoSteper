import os

from autosteper import AutoSteper


para = {
    'pristine_path': r'geom.xyz',
    'root': r'./',
    'gen_para': {'group': 'CF3',
                 'geom_mode': 'pre_defined',
                 'gen_core_path': r'I:\new_nauty\my_v\usenauty\bin\cagesearch.exe'
                 },
    'opt_mode': 'xtb',
    'opt_para': {
        'cmd_list': [r'/home/mkliu/anaconda3/envs/env001/bin/xtb', '--opt', 'tight', '--json', '--cycles 10'],
        'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
        'deal_wrong_mode': 'Report',
        'mach_para': {
            'batch_type': "Torque",
            'context_type': "SSHContext",
            'remote_root': '/home/mkliu/test_dpdispatcher/',
            'remote_profile': {
                "hostname": "219.245.39.76",
                "username": "mkliu",
                "password": "mkliu123",
                "port": 22,
                "timeout": 10
            }
        },
        'resrc_para': {
            'number_node': 1,
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
            'sub_batch_size': 48
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
