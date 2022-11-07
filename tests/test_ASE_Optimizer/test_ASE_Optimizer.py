from autosteper.optimizers import switch_optimizers

para = {
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
    'ase_para': {'model_path': r'last.ckpt',
                 'py_script': r'parallel_unit.py',
                 'num_pll': 5,
                 'base_node': 0,
                 'cpu_per_worker': 6},
}

an_optimizer = switch_optimizers(mode='ase', para=para)
# raw folder is pre-defined
an_optimizer.set_folders()
opt_status = an_optimizer.opt()
assert opt_status == 0
