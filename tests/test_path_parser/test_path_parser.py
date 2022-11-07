import os

from autosteper import AutoSteper
from autosteper.cage import Cage
from autosteper.parser import Path_Parser

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

a_cage = Cage(pristine_path=r'./geom.xyz')
a_cage.set_workbase(root=r'/home/mkliu/schnet_opt/paper_4_27/10_6')

path_para = {
    'step': 1,
    'start': 1,
    'q_add_num': 4,
    'q_path_rank': 10,
    'q_isomer_rank': 1,
    'log_low_e_num': 10,

}

a_path_parser = Path_Parser(path_para=path_para,
                            dump_root=r'/home/mkliu/schnet_opt/paper_4_27/10_6/parser',
                            q_cage=a_cage)
a_path_parser.get_path_info()