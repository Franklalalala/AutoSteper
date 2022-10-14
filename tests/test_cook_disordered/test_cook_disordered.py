import os

from autosteper.parser import cook_disordered


cook_disordered(disordered_root=r'./disordered_logs',
                dump_root=r'./',
                keep_top_k_pathway=500,
                step=1,
                log_mode='gauss')
