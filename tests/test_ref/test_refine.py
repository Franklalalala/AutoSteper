import os

from autosteper.parser import refine

mach_para = {
    'batch_type': "Torque",
    'context_type': "LocalContext",
    'remote_root': '/home/mkliu/test_dpdispatcher/',
    'remote_profile': {
        None
    }
}

resrc_para = {
    'number_node': 1,
    'cpu_per_node': 6,
    'gpu_per_node': 0,
    'group_size': 25,
    'queue_name': "batch",
    'envs': {
        "OMP_STACKSIZE": "4G",
        "OMP_NUM_THREADS": "3,1",
        "OMP_MAX_ACTIVE_LEVELS": "1",
        "MKL_NUM_THREADS": "3"
    },
    'sub_batch_size': 200
}

para = {
    'start': 1,
    'stop': 12,
    'step': 1,
    'opt_mode': 'xtb',
    'opt_para': {
        'cmd_list': [r'/home/mkliu/anaconda3/envs/env001/bin/xtb', '--opt', 'tight', '--json'],
        'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
        'deal_wrong_mode': 'Report',
        'mach_para': mach_para,
        'resrc_para': resrc_para,
    },
    'cutoff': {
        'mode': 'value_and_rank',
        'value': 1,
        'rank': 200
    }
}

refine(old_workbase=r'/home/mkliu/schnet_opt/paper_4_27/10_6/geom',
       new_workbase=r'/home/mkliu/schnet_opt/paper_4_27/10_6/test_ref/geom',
       ref_para=para)
