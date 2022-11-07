from autosteper import XTB_Optimizer


an_optimizer = XTB_Optimizer(opt_para={
    'cmd_list': [r'/home/mkliu/anaconda3/envs/env001/bin/xtb', '--opt', 'tight', '--json'],
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
        'number_node': 6,
        'cpu_per_node': 6,
        'gpu_per_node': 0,
        'group_size': 10,
        'queue_name': "batch",
        'sub_batch_size': 50
    },
})

# raw folder is pre-defined
an_optimizer.set_folders()
opt_status = an_optimizer.opt()
assert opt_status == 0
