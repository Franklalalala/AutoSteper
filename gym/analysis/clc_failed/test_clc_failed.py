import os

from autosteper.parser import clc_failed


failed_counts, add_list = clc_failed(workbase=r'dummy_Br_workbase',
                                     dump_pic_path=r'failed_distribution.png')

with open('swr_info.txt', 'w') as f_w:
    f_w.write('Addon list: ')
    f_w.write(str(add_list))
    f_w.write('\n========================================\n')
    for a_key in failed_counts:
        f_w.write(f'Counts:     {str(failed_counts[a_key])},   Failed type: {str(a_key)}')
        f_w.write('\n========================================\n')