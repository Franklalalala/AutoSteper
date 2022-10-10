import pandas as pd
from ase.io import read
from autosteper.parser import simple_test_iso, simple_log_relatives, strict_scatter_relatives


target_info_path = r'/home/mkliu/schnet_opt/paper_4_27/10_6/test_ref/geom/10addons/passed_info.pickle'
atoms = read('C60Cl10_exp.xyz')
rank = simple_test_iso(q_atoms=atoms, passed_info_path=target_info_path, top_k=10)
target_xyz_path = pd.read_pickle(target_info_path)['xyz_path'][rank - 1]
print(target_xyz_path)
q_atoms = read(target_xyz_path)

simple_log_relatives(workbase=r'/home/mkliu/schnet_opt/paper_4_27/10_6/test_ref/geom',
                     dump_log_path=r'/home/mkliu/schnet_opt/paper_4_27/10_6/test_ref/rel.log',
                     step=1,
                     fst_add_num=1,
                     final_add_num=12,
                     q_atoms=q_atoms)

strict_scatter_relatives(workbase=r'/home/mkliu/schnet_opt/paper_4_27/10_6/test_ref/geom',
                         dump_folder=r'/home/mkliu/schnet_opt/paper_4_27/10_6/test_ref/dump',
                         step=1,
                         fst_add_num=1,
                         final_add_num=12,
                         q_atoms=q_atoms)
