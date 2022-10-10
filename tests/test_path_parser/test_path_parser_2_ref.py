import os

from autosteper import AutoSteper
from autosteper.cage import Cage
from autosteper.parser import Path_Parser

a_cage = Cage(pristine_path=r'./geom.xyz')
a_cage.set_workbase(root=r'/home/mkliu/schnet_opt/paper_4_27/10_6')

mach_para = {
    'batch_type': "Torque",
    'context_type': "LocalContext",
    'remote_root': '/home/mkliu/test_dpdispatcher/',
    'remote_profile': {
        None
    }
}

resrc_para = {
    'number_node': 2,
    'cpu_per_node': 6,
    'gpu_per_node': 0,
    'group_size': 5,
    'queue_name': "batch",
    'envs': {
        "OMP_STACKSIZE": "4G",
        "OMP_NUM_THREADS": "3,1",
        "OMP_MAX_ACTIVE_LEVELS": "1",
        "MKL_NUM_THREADS": "3"
    },
    'sub_batch_size': 40
}

refiner_para = {
    'has_parity': True,
    'cage_size': a_cage.size,
    'refine_top_k': 15,
    'opt_mode': 'xtb',
    'cmd_list': [r'/home/mkliu/anaconda3/envs/env001/bin/xtb', '--opt', 'tight', '--json'],
    'out_list': ['xtbopt.xyz', 'xtbopt.log', 'xtbout.json'],
    'deal_wrong_mode': 'Report',
    'mach_para': mach_para,
    'resrc_para': resrc_para
}

path_para = {
    'is_mix': True,
    'step': 1,
    'start': 1,
    'q_add_num': 4,
    'q_path_rank': 20,
    'q_isomer_rank': 5,
    'log_low_e_num': 10,

}

a_path_parser = Path_Parser(path_para=path_para,
                            dump_root=r'/home/mkliu/schnet_opt/paper_4_27/10_6/parser2_ref',
                            q_cage=a_cage,
                            refiner_para=refiner_para)
a_path_parser.get_path_info()
