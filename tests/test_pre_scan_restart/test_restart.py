import os

from autosteper import AutoSteper
from somenet import net
from somenet import calculator
import torch


model_path = r'xx/Cl_final.ckpt'
state_dict = torch.load(model_path)
net.load_state_dict(state_dict["state_dict"])
calculator = calculator(net=net)

para = {
    'pristine_path': r'geom.xyz',
    'root': r'./',
    'gen_para': {'group': 'Cl',
                 'geom_mode': 'pre_defined',
                 'gen_core_path':  r"/home/mkliu/nauty/usenauty/bin/cagesearch",
                 },
    'opt_mode': 'ase',
    'opt_para': {
        'cmd_list': [r'/home/mkliu/anaconda3/envs/molnet/bin/python3.8'],
        'out_list': ['cooked'],
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
        # A parallel distribution file is needed.
        'ase_para': {'model_path': r'xx/last.ckpt',
                     'py_script': r'xx/parallel_unit.py',
                     'num_pll': 8,
                     'base_node': 0,
                     'cpu_per_worker': 6},
    },
    'run_para': {
        'start': 1,
        'stop': 6,
        'step': 1,
        'wht_list_para': {
            'mode': 'rank',
            'rank': 10
        }
    },
    'pre_scan_para': {
        'ps_num_list': [2, 3],
        'calculator': calculator,
        'ps_cut_para': {
            'mode': 'rank',
            'rank': 80
        }
    }
}

auto = AutoSteper(para=para)
auto.restart(restart_add_num=5)
