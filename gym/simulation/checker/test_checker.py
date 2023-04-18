import os

from ase.io import read
from autosteper.cage import Cage
from autosteper.checker import Checker


pristine_cage = Cage(pristine_path=r'C60_000001812opted.xyz')
a_checker = Checker(pst_cage=pristine_cage, group='Cl')
os.chdir('error_xyz')
for i in range(1, 7):
    os.chdir(str(i))
    final_image = read('opt.xyz')
    for a_file in os.listdir('./'):
        if a_file.endswith('.xyz') and a_file != 'opt.xyz':
            q_name = os.path.splitext(a_file)[0]
    job_status = a_checker.check_a_job(q_atoms=final_image, q_name=q_name)
    assert job_status == i
    os.chdir('./..')

os.chdir('./7')
pristine_cage = Cage(pristine_path=r'pristine_cage.xyz')
a_checker = Checker(pst_cage=pristine_cage, group='CF3')
final_image = read('xtbopt.xyz')
q_name = r'0000004zu6td'
job_status = a_checker.check_a_job(q_atoms=final_image, q_name=q_name)
assert job_status == 7
