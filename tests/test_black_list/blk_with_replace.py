import os
import time

import torch
from ase.optimize import *
from autosteper.Autosteper import AutoSteper
from somenet import calculator
from somenet import net

model_path = r'xx/Cl_final.ckpt'
state_dict = torch.load(model_path)
net.load_state_dict(state_dict["state_dict"])
calculator = calculator(net=net)

para = {
    'workbase': r'xx',
    'pristine_cage': r'xx/C106.xyz',
    'gen_para': {
        'gen_core_path': r"xx/cagesearch",
        'geom_mode': 'pre_defined',
        'group': 'Cl',
        'skin': 0.15
    },
    'opt_mode': 'ase',
    'opt_para': {
        'calculator': calculator,
        'is_Opt_Twice': False,
        'init_cycle': 100,
        'fmax': 0.0005,
        'deal_wrong_mode': 'report',
        'ase_optimizer': FIRE,
        'is_pll': True,
        'pll_para': {
            'num_worker': 9,
            'cpu_per_worker': 5,
            'base_node': 5,
            'pll_unit_file': r'xx/parallel_unit.py'
        },
    },
    'run_para': {
        'start': 1,
        'step': 1,
        'stop': 24,
        'wht_list_para': {
            'mode': 'value_rank',
            'value': 1,
            'rank': 200
        }
    },
    'blk_para': {
        'start_clct_num': 3,
        'end_chk_num': 21,
        'clct_unstb': True,
        'unstb_para': {
            'container_size': 2,
            'mode': 'value_rank',
            'value': 0.3,
            'rank': 500
        }
    },
}

auto = AutoSteper(para)
auto.run()

para = {
    'workbase': r'xx',
    'pristine_cage': r'xx/C106.xyz',
    'gen_para': {
        'gen_core_path': r"xx/cagesearch",
        'geom_mode': 'pre_defined',
        'group': 'Cl',
        'skin': 0.15
    },
    'opt_mode': 'ase',
    'opt_para': {
        'calculator': calculator,
        'is_Opt_Twice': False,
        'init_cycle': 100,
        'fmax': 0.0005,
        'deal_wrong_mode': 'report',
        'ase_optimizer': FIRE,
        'is_pll': True,
        'pll_para': {
            'num_worker': 9,
            'cpu_per_worker': 5,
            'base_node': 5,
            'pll_unit_file': r'xx/parallel_unit.py'
        },
    },
    'run_para': {
        'start': 1,
        'step': 1,
        'stop': 24,
        'wht_list_para': {
            'mode': 'value_rank',
            'value': 1,
            'rank': 200
        }
    },
    'blk_para': {
        'start_clct_num': 3,
        'end_chk_num': 21,
        'clct_unstb': True,
        'unstb_para': {
            'container_size': 1,
            'mode': 'value_rank',
            'value': 0.3,
            'rank': 500
        }
    },
}

auto = AutoSteper(para)
auto.run()
