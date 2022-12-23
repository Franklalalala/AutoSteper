import os
from autosteper.Autosteper import AutoSteper
from ase.optimize import *
from torchlightmolnet.caculator import torchCaculator
from torchlightmolnet.lightning.molnet import LightMolNet
from torchlightmolnet.dataset.atomref import refatoms_xTB, get_refatoms
from torchlightmolnet import Properties
import torch
import time


para = {
    'workbase': r'/home/mkliu/dummy_test_floder/C66H4/workbase',
    'pristine_cage': r'/home/mkliu/dummy_test_floder/C66H4/xyz/C66_000004466opt.xyz',
    'gen_para': {
        'gen_core_path': r"/home/mkliu/nauty/usenauty/bin/cagesearch",
        'geom_mode': 'pre_defined',
        'group': 'H',
        'skin': 0.15
    },
    'opt_mode': 'xtb',
    'opt_para': {
        'xtb_path':  r'/home/mkliu/anaconda3/envs/env001/bin/xtb',
        'cmd_list': ['--opt', 'tight', '--json'],
        'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
        'is_Opt_Twice': False,
        'init_cycle': 100,
        'deal_wrong_mode': 'report',
        'mach_para': {
            'batch_type': "Torque",
            'context_type': "LocalContext",
            'remote_root': '/home/mkliu/test_dpdispatcher/',
            'remote_profile': {
                None
            }
        },
        'resrc_para': {
            'number_node': 1,
            'cpu_per_node': 6,
            'gpu_per_node': 0,
            'group_size': 125,
            'queue_name': "batch",
            'envs': {
                "OMP_STACKSIZE": "4G",
                "OMP_NUM_THREADS": "3,1",
                "OMP_MAX_ACTIVE_LEVELS": "1",
                "MKL_NUM_THREADS": "3"
            }
        },
        'sub_batch_size': 1000
    },
    'run_para': {
        'start': 1,
        'step': 1,
        'stop': 12,
        'wht_list_para': {
            'mode': 'value_rank',
            'value': 1,
            'rank': 200
        }
    },
    'path_para': {
        'q_path_num': 10,
        'q_add_num': 12,
        'q_low_e_num': 1,
        'log_low_e_num': 100,
        # 'ctl_path_para': {
        #     'max_path_num' : 10000,
        #     'ctl_parent_num': 2
        # }
    }
}


auto = AutoSteper(para)
auto.run()
auto.path_parser.get_path_info()
