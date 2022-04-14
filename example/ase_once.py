import os
from autosteper.Autosteper import AutoSteper
from ase.calculators.emt import EMT
from ase.optimize import *




para = {
    'workbase': r'F:\AutoSteper\example',
    'pristine_cage': r'F:\AutoSteper\example\c50.xyz',

    'gen_para': {
        'gen_core_path': r"I:\new_nauty\my_v\usenauty\bin\cagesearch.exe",
        'geom_mode': 'pre_defined',
        'group': 'H',
        'skin': 0.15
    },
    'opt_mode': 'ase',
    'opt_para': {
        'calculator': EMT(),
        'is_Opt_Twice': False,
        'init_cycle': 100,
        'fmax': 0.0005,
        'deal_wrong_mode': 'report',
        'ase_optimizer': FIRE
    },
    'run_para': {
        'start': 1,
        'step': 1,
        'stop': 2,
        'cutoff_mode': 'value_rank',
        'cutoff_value': 0.037,
        'cutoff_rank': 1
    },
    'path_para': {
        'q_path_num': 200,
        'q_add_num': 2,
        'q_low_e_num': 65,
        'log_low_e_num': 3,
    }
}


auto = AutoSteper(para)
auto.run()
auto.path_parser.get_path_info()
