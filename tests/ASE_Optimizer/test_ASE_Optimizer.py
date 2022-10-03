from autosteper.optimizers import switch_optimizers

para = {
    'cmd_list': [r'/home/mkliu/anaconda3/envs/molnet/bin/python3.8'],
    'out_list': ['cooked'],
    'deal_wrong_mode': 'Report',
    'mach_para': {
        'batch_type': "Torque",
        'context_type': "SSHContext",
        'remote_root': '/home/mkliu/test_dpdispatcher/',
        'remote_profile': {
            "hostname": "219.245.39.76",
            "username": "mkliu",
            "password": "mkliu123",
            "port": 22,
            "timeout": 10
        }
    },
    'resrc_para': {
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
        'sub_batch_size': 5
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
