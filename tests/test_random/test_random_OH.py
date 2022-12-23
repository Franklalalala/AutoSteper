import os

from autosteper import AutoSteper


para = {
    'pristine_path': r'/home/mkliu/schnet_opt/12_9_OH/data/geom.xyz',
    'root': r'/home/mkliu/schnet_opt/12_9_OH/data/dump',
    'gen_para': {'group': 'OH',
                 'geom_mode': 'pre_defined',
                 'gen_core_path': r'/home/mkliu/nauty/usenauty/bin/cagesearch'
                 },
    'opt_mode': 'xtb',
    'opt_para': {
        'has_parity': True,
        'cage_size': 60,
        'cmd_list': [r'/home/mkliu/anaconda3/envs/env001/bin/xtb', '--opt', '--json', '--cycles 50'],
        'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
        'deal_wrong_mode': 'Report',
        'mach_para': {
            'batch_type': "Torque",
            'context_type': "LocalContext",
            'remote_root': '/home/mkliu/test_dpdispatcher/',
            'remote_profile': {
                None
            }
        },
        'resrc_para': {
            'number_node': 3,
            'cpu_per_node': 6,
            'gpu_per_node': 0,
            'group_size': 75,
            'queue_name': "batch",
            'envs': {
                "OMP_STACKSIZE": "4G",
                "OMP_NUM_THREADS": "3,1",
                "OMP_MAX_ACTIVE_LEVELS": "1",
                "MKL_NUM_THREADS": "3"
            }
        },
    },
    'random_para': {
        'addon_list': list(range(4, 53, 3)),
        'random_num': 300,
        'try_times': 3
    }
}

# raw folder is pre-defined
auto = AutoSteper(para=para)
auto.random()
