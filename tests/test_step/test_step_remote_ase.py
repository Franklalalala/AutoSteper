import os

from autosteper import AutoSteper


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
        'group_size': 1,
        'queue_name': "batch",
        'envs': {
            "OMP_STACKSIZE": "4G",
            "OMP_NUM_THREADS": "3,1",
            "OMP_MAX_ACTIVE_LEVELS": "1",
            "MKL_NUM_THREADS": "3"
        },
        'sub_batch_size': 8
    }



para = {
    'pristine_path': r'./geom.xyz',
    'root': r'/home/mkliu/schnet_opt/paper_4_27/10_6',
    'gen_para': {'group': 'Cl',
                 'geom_mode': 'pre_defined',
                 'gen_core_path':  r"/home/mkliu/nauty/usenauty/bin/cagesearch",
                 },
    'opt_mode': 'ase',
    'opt_para': {
        'cmd_list': [r'/home/mkliu/anaconda3/envs/molnet/bin/python3.8'],
        'out_list': ['cooked'],
        'deal_wrong_mode': 'Report',
        'mach_para': mach_para,
        'resrc_para': resrc_para,
        'ase_para': {'model_path': r'/home/mkliu/schnet_opt/paper_4_27/10_6/last.ckpt',
                     'py_script': r'/home/mkliu/schnet_opt/paper_4_27/10_6/parallel_unit.py',
                     'num_pll': 8,
                     'base_node': 0,
                     'cpu_per_worker': 6},
    },
    'run_para': {
        'start': 1,
        'stop': 12,
        'step': 1,
        'wht_list_para': {
            'mode': 'value_and_rank',
            'value': 1,
            'rank': 200
        }
    },
}

# raw folder is pre-defined
auto = AutoSteper(para=para)
auto.run()
