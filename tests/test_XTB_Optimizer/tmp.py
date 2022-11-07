import pandas as pd
import shutil


# info = pd.read_pickle(r'F:\AutoSteper_new\AutoSteper\tests\test_random\geom_v2\55addons\status_info.pickle')
# info = pd.read_pickle(r'F:\AutoSteper_new\AutoSteper\tests\test_step\geom\1addons\passed_info.pickle')
# info = pd.read_pickle('status_info.pickle')
info = pd.read_pickle('passed_info.pickle')
print(info['name'])
print(info['energy'])

# info_2 = pd.read_pickle(r'F:\AutoSteper\tests\test_step\geom\1addons\parent_info.pickle')
# shutil.copytree(dst=r'F:\AutoSteper_new\AutoSteper\tests\XTB_Optimizer\cooking\H2_0',
#                 src=r'F:\AutoSteper_new\AutoSteper\tests\deal_wrong_recursive\cooking\H2_0',
#                 dirs_exist_ok=True)

# info.index = list(range(len(info.index)))
a=1