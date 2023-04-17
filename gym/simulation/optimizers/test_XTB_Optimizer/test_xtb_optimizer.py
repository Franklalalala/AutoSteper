from autosteper import XTB_Optimizer


an_optimizer = XTB_Optimizer(opt_para={
    'cmd_list': [r'path/to/xtb', '--opt', 'tight', '--json'],
    'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
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
        'number_node': 1,
        'cpu_per_node': 6,
        'gpu_per_node': 0,
        'group_size': 2,
        'queue_name': "xx",
        'sub_batch_size': 10,
        'envs': {
            'dummy': 'dummy/path',
            # source your env here
            "OMP_STACKSIZE": "4G",
            "OMP_NUM_THREADS": "3,1",
            "OMP_MAX_ACTIVE_LEVELS": "1",
            "MKL_NUM_THREADS": "3"
        }
    },
})

# raw folder is pre-defined
an_optimizer.set_folders()
opt_status = an_optimizer.opt()
assert opt_status == 0
