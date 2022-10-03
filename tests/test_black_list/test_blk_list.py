import os

import torch
from autosteper import AutoSteper
from torchlightmolnet import Properties
from torchlightmolnet.caculator import torchCaculator
from torchlightmolnet.dataset.atomref import refatoms_xTB, get_refatoms
from torchlightmolnet.lightning.molnet import LightMolNet

model_path = r'/home/mkliu/schnet_opt/paper_4_27/Cl_final.ckpt'
net = LightMolNet(atomref=get_refatoms(refatoms_xTB)[Properties.energy_U0])
state_dict = torch.load(model_path)
net.load_state_dict(state_dict["state_dict"])
calculator = torchCaculator(net=net)

para = {
    'pristine_path': r'C108.xyz',
    'root': r'./',
    'gen_para': {'group': 'Cl',
                 'geom_mode': 'pre_defined',
                 'gen_core_path': r"/home/mkliu/nauty/usenauty/bin/cagesearch",
                 },
    'opt_mode': 'ase',
    'opt_para': {
        'cmd_list': [r'/home/mkliu/anaconda3/envs/molnet/bin/python3.8'],
        'out_list': ['cooked'],
        'deal_wrong_mode': 'Report',
        'mach_para': {
            'batch_type': "Torque",
            'context_type': "LocalContext",
            'remote_root': '/home/mkliu/test_dpdispatcher/',
            'remote_profile': None
        },
        'resrc_para': {
            'number_node': 1,
            'cpu_per_node': 6,
            'gpu_per_node': 0,
            'group_size': 1,
            'queue_name': "batch",
            'envs': {},
            'sub_batch_size': 8
        },
        # A parallel distribution file is needed.
        'ase_para': {'model_path': r'/home/mkliu/anaconda3/envs/molnet/AutoSteper/tests/test_pre_scan/last.ckpt',
                     'py_script': r'/home/mkliu/anaconda3/envs/molnet/AutoSteper/tests/test_pre_scan/parallel_unit.py',
                     'num_pll': 8,
                     'base_node': 0,
                     'cpu_per_worker': 6},
    },
    'run_para': {
        'start': 1,
        'stop': 10,
        'step': 1,
        'wht_list_para': {
            'mode': 'rank',
            'rank': 20
        }
    },
    'blk_para': {
        'start_clct_num': 2,
        'final_chk_num': 8,
        'clct_unstb': True,
        'unstb_para': {
            'mode': 'rank',
            'rank': 10,
        },
        'container_size': 3
    },
    'pre_scan_para': {
        'start_ps_para': 2,
        'final_ps_para': 8,
        'calculator': calculator,
        'ps_cut_para': {
            'mode': 'rank',
            'rank': 200
        }
    }
}

auto = AutoSteper(para=para)
auto.run()
