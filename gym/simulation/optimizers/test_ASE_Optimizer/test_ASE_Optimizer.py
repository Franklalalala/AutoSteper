from autosteper.optimizers import switch_optimizers

para = {
    'cmd_list': [r'path/to/python3.8'],
    'out_list': ['cooked'],
    'deal_wrong_mode': 'Report',
    'mach_para': {
        'batch_type': "Torque",
        'context_type': "LocalContext",
        'remote_root': 'path/to/test_dpdispatcher',
        'remote_profile': None,
    },
    'resrc_para': {
        'number_node': 6,
        'cpu_per_node': 6,
        'gpu_per_node': 0,
        'group_size': 1, # This parameter needs to stay with 1 for ase optimizers
        'queue_name': "batch",
        'envs': {
            "OMP_STACKSIZE": "4G",
            "OMP_NUM_THREADS": "3,1",
            "OMP_MAX_ACTIVE_LEVELS": "1",
            "MKL_NUM_THREADS": "3"
        }
    },
    # A parallel distribution file is needed.
    'ase_para': {'model_path': r'path/to/Cl_final.ckpt',
                 'py_script': r'path/to/parallel_unit.py',
                 'num_pll': 8,
                 'base_node': 0,
                 'cpu_per_worker': 3},
}

an_optimizer = switch_optimizers(mode='ase', para=para)
# raw folder is pre-defined
an_optimizer.set_folders()
opt_status = an_optimizer.opt()
assert opt_status == 0
