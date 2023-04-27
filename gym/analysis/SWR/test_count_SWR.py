import os

from autosteper.parser import count_SWR


swr_clctor, add_list = count_SWR(swr_1_workbase=r'./11_12/q_11_to_tgt_12',
                                 swr_2_workbase=r'./11_12/q_12_to_tgt_11',
                                 swr_1_legend=r'$\rm ^{\#11}C_{84}$ to $\rm ^{\#12}C_{84}$',
                                 swr_2_legend=r'$\rm ^{\#12}C_{84}$ to $\rm ^{\#11}C_{84}$',
                                 dump_pic_path=r'./swr_count_result.png')

with open('swr_info.txt', 'w') as f_w:
    f_w.write('Addon list:\n')
    f_w.write(str(add_list))
    f_w.write('\n========================================\n')
    for a_key in swr_clctor:
        f_w.write('Name:\n')
        f_w.write(str(a_key))
        f_w.write('\n')
        f_w.write('Counts:\n')
        f_w.write(str(swr_clctor[a_key]))
        f_w.write('\n========================================\n')
