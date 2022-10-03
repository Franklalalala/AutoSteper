from autosteper import Multi_Optimizer


mach_para = {
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
    }

resrc_para = {
        'number_node': 6,
        'cpu_per_node': 6,
        'gpu_per_node': 0,
        'group_size': 2,
        'queue_name': "batch",
        'envs': {
            "OMP_STACKSIZE": "4G",
            "OMP_NUM_THREADS": "3,1",
            "OMP_MAX_ACTIVE_LEVELS": "1",
            "MKL_NUM_THREADS": "3"
        },
        'sub_batch_size': 10
    }

xtb_para_1 = {
    'opt_mode': 'xtb',
    'cmd_list': [r'/home/mkliu/anaconda3/envs/env001/bin/xtb', '--opt', 'tight', '--json', '--cycles 15'],
    'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
    'deal_wrong_mode': 'Report',
    'mach_para': mach_para,
    'resrc_para': resrc_para,
}

gau_para = {
    'opt_mode': 'gaussian',
    'cmd_list': [r'g16'],
    'out_list': ['gau.log', 'gau.chk'],
    'deal_wrong_mode': 'Report',
    'mach_para': mach_para,
    'resrc_para': resrc_para,
    'gau_para': {
        'cmd_line': r'opt b3lyp/3-21G',
        'charge': 0,
        'multi': 1,
        'cpu_per_worker': 6,
        'mem': '10GB'
    },
    'wht_list_para': {
        'mode': 'rank',
        'rank': 5
    }
}

xtb_para_2 = {
    'opt_mode': 'xtb',
    'cmd_list': [r'/home/mkliu/anaconda3/envs/env001/bin/xtb', '--opt', 'tight', '--json'],
    'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
    'deal_wrong_mode': 'Report',
    'mach_para': mach_para,
    'resrc_para': resrc_para,
    'wht_list_para': {
        'mode': 'value',
        'value': 0.1,
        'nimg_th': 2
    }
}

# raw folder is pre-defined
an_optimizer = Multi_Optimizer(opt_para=[xtb_para_1, gau_para, xtb_para_2])
opt_status = an_optimizer.opt()
assert opt_status == 0
