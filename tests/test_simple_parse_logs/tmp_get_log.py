import os
import shutil


dump = r'/home/mkliu/dummy_test_floder/C66H4/workbase/temp/log'
old_old_root = r'/home/mkliu/dummy_test_floder/C66H4/workbase/C66_000004169opt'
for ii in [3,6,8]:
    old_root = os.path.join(old_old_root, f'{ii}addons', 'opt_100')
    for i in os.listdir(old_root)[:5]:
        os.chdir(old_root)
        os.chdir(i)
        shutil.copy(src='xtbopt.log', dst=os.path.join(dump, f'{i}.log'))

