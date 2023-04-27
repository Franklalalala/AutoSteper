# The cook_disordered_results is cook_disordered function results, which is highly structured and informative.
# proj_ring_seq is defined to make a better visualization results.
# One may try_every_hexagon in plotter section to pick a ring in advance.

import os

from autosteper.parser import swr_pipeline


swr_pipeline(sorted_root_1=r'./cook_disordered_results/11/sorted',
             sorted_root_2=r'./cook_disordered_results/12/sorted',
             is_unique=True,
             is_low_e=True,
             label_1='11',
             label_2='12',
             legend_swr_1_to_2=r'$\rm ^{\#11}C_{84}$ to $\rm ^{\#12}C_{84}$',
             legend_swr_2_to_1=r'$\rm ^{\#12}C_{84}$ to $\rm ^{\#11}C_{84}$',
             dump_root=r'./11_12',
             proj_ring_seq_1=[71, 72, 73, 74, 75, 76],
             proj_ring_seq_2=[71, 72, 73, 74, 75, 76],
             group_symbol='Cl',
             dpi=100)
