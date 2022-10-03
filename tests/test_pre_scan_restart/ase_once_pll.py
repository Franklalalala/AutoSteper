import os
from autosteper.Autosteper import AutoSteper
from ase.optimize import *
from torchlightmolnet.caculator import torchCaculator
from torchlightmolnet.lightning.molnet import LightMolNet
from torchlightmolnet.dataset.atomref import refatoms_xTB, get_refatoms
from torchlightmolnet import Properties
import torch


model_path = r'/home/mkliu/schnet_opt/paper_4_27/Cl_final.ckpt'
net = LightMolNet(atomref=get_refatoms(refatoms_xTB)[Properties.energy_U0])
state_dict = torch.load(model_path)
net.load_state_dict(state_dict["state_dict"])
calculator = torchCaculator(net=net)


para = {
    'workbase': r'/home/mkliu/schnet_opt/test_pre_scan',
    'pristine_cage': r'/home/mkliu/schnet_opt/paper_4_27/xyz/C50_271_saturn.xyz',
    'gen_para': {
        'gen_core_path': r"/home/mkliu/nauty/usenauty/bin/cagesearch",
        'geom_mode': 'pre_defined',
        'group': 'Cl',
        'skin': 0.15
    },
    'is_pre_scan': True,
    'pre_scan_para': {
        'calculator': calculator,
        'start_ps_num': 2,
        'ps_cut_para':{
            'mode': 'rank',
            'rank': 50
        },
    },
    'opt_mode': 'ase',
    'opt_para': {
        'calculator': calculator,
        'is_Opt_Twice': False,
        'init_cycle': 100,
        'fmax': 0.0005,
        'deal_wrong_mode': 'report',
        'ase_optimizer': FIRE,
        # 'is_pll': False
        'is_pll': True,
        'pll_para': {
            'num_worker': 11,
            'cpu_per_worker': 6,
            'base_node': 1,
            'pll_unit_file': r'/home/mkliu/schnet_opt/paper_4_27/parallel/parallel_unit.py'
        },
    },
    'run_para': {
        'start': 2,
        'step': 2,
        'stop': None,
        'cutoff_mode': 'value_rank',
        'cutoff_value': 0.037,
        'cutoff_rank': 200
    },
    'path_para': {
        'q_path_num': 200,
        'q_add_num': 12,
        'q_low_e_num': 65,
        'log_low_e_num': 3,
    }
}


auto = AutoSteper(para)
auto.run()
auto.path_parser.get_path_info()
