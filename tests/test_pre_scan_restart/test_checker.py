from autosteper.checker import Checker
from autosteper.cage import Cage
from ase.io import read

# a_cage = Cage(pristine_path=r'/home/mkliu/anaconda3/envs/molnet/AutoSteper/tests/test_pre_scan/geom.xyz')
# a_checker = Checker(pst_cage=a_cage, group='Cl')
# a_checker.check(passed_info_path=r'passed_info.pickle', status_info_path=r'status_info.pickle')


a_cage = Cage(pristine_path=r'F:\AutoSteper_new\AutoSteper\tests\test_random\geom.xyz')
a_checker = Checker(pst_cage=a_cage, group='Cl')

# a_checker.check(passed_info_path=r'passed_info.pickle', status_info_path=r'status_info.pickle')


q_atoms = read(r'xtbopt.xyz')
a_status = a_checker.check_a_job(q_atoms=q_atoms, q_name='8ckk6x9pwq0z')
print(a_status)