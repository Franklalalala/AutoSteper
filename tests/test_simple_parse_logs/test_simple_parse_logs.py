import os

from autosteper.tools import simple_parse_logs

simple_parse_logs(dump_root=r'F:\AutoSteper\tests\test_simple_parse_logs\gauss',
                  src_root=r'F:\AutoSteper\tests\test_simple_parse_logs\gauss\disordered_logs',
                  group='OH',
                  cage_size=60,
                  mode='gauss')

simple_parse_logs(dump_root=r'F:\AutoSteper\tests\test_simple_parse_logs\xtb',
                  src_root=r'F:\AutoSteper\tests\test_simple_parse_logs\xtb\disordered_logs',
                  mode='xyz')
