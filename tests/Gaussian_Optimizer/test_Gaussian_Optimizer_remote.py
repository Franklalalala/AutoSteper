from autosteper import Gaussian_Optimizer

an_optimizer = Gaussian_Optimizer(opt_para={
    'cmd_list': [r'g16'],
    'out_list': ['gau.log', 'gau.chk'],
    'deal_wrong_mode': 'Report',
    'mach_para': {
        'batch_type': "Torque",
        'context_type': "LocalContext",
        'remote_root': '/home/mkliu/test_dpdispatcher/',
        'remote_profile': None
    },
    'resrc_para': {
        'number_node': 1,
        'cpu_per_node': 9,
        'gpu_per_node': 0,
        'group_size': 2,
        'queue_name': "batch",
        'envs': {
            "OMP_STACKSIZE": "4G",
            "OMP_NUM_THREADS": "3,1",
            "OMP_MAX_ACTIVE_LEVELS": "1",
            "MKL_NUM_THREADS": "3"
        },
        'sub_batch_size': 25
    },
    'gau_para': {
        'cmd_line': r'opt b3lyp/3-21G',
        'charge': 0,
        'multi': 1,
        'cpu_per_worker': 6,
        'mem': '10GB'
    },
})

# raw folder is pre-defined
an_optimizer.set_folders()
opt_status = an_optimizer.opt()
assert opt_status == 0
