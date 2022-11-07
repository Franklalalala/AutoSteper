import os

from autosteper.tools import simple_parse_logs

simple_parse_logs(dump_root=r'./gauss',
                  src_root=r'./gauss/disordered_logs',
                  group='OH',
                  cage_size=60,
                  mode='gauss')

simple_parse_logs(dump_root=r'./xtb',
                  src_root=r'./xtb/disordered_logs',
                  mode='xyz')
