from autosteper import Gaussian_Optimizer

an_optimizer = Gaussian_Optimizer(opt_para={
    'cmd_list': [r'g16'],
    'out_list': ['gau.log', 'gau.chk'],
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
