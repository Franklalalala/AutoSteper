import os
import shutil
import sys
import warnings
import networkx as nx
from pathlib import Path
from typing import Optional, Union
from typing import Union
import math
import dpdata
import numpy as np
import pandas as pd
from ase import Atoms, Atom
from ase.io import read, write
from ase.neighborlist import build_neighbor_list
from ase.visualize import view
from fullerenedatapraser.io.recursion import recursion_files
from fullerenedatapraser.io.xyz import simple_read_xyz_xtb
from fullerenedatapraser.molecular.fullerene import FullereneFamily
from tqdm import tqdm
from ase.optimize import BFGS, LBFGS
from ase.optimize.optimize import Optimizer
from ase.calculators.emt import EMT
from ase.io import Trajectory
from torchlightmolnet.caculator import torchCaculator
from torchlightmolnet.lightning.molnet import LightMolNet
from torchlightmolnet.dataset.atomref import refatoms_xTB, get_refatoms
from torchlightmolnet import Properties
import torch
import subprocess
from ase.units import Hartree, eV
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from dpdispatcher import Task, Submission, Machine, Resources


def lazy_rm(a_path):
    if os.path.exists(a_path):
        shutil.rmtree(a_path)
    os.makedirs(exist_ok=True, name=a_path)
    return os.path.abspath(a_path)

dump = lazy_rm(r'./dump')
dst = lazy_rm(r'./opt')
xyz_path = os.path.abspath(r'./xyz')


def run_a_batch(path_source: str, path_destination: str, cmd_list: list, out_list: list, sub_batch_size: int=None):
    path_destination = os.path.abspath(path_destination)
    mach = Machine(batch_type="Torque",
                        context_type="LocalContext",
                        remote_root=r'/home/mkliu/test_dpdispatcher/',
                        remote_profile=None,
                        local_root=path_destination)
    resrc = Resources(number_node=1,
                           cpu_per_node=8,
                           gpu_per_node=0,
                           group_size=6,
                           queue_name='batch',
                           envs={
                "OMP_STACKSIZE": "4G",
                "OMP_NUM_THREADS": "3,1",
                "OMP_MAX_ACTIVE_LEVELS": "1",
                "MKL_NUM_THREADS": "3"
            }
                           )
    task_list = []
    for item in os.listdir(path_source):
        isomer_name = os.path.splitext(item)[0]
        os.makedirs(os.path.join(path_destination, isomer_name), exist_ok=True)
        path_item = os.path.join(path_source, item)
        shutil.copy(path_item, os.path.join(path_destination, isomer_name, item))
        cmd_list.append(item)
        a_task = Task(command=' '.join(cmd_list),
                      task_work_path=f"{isomer_name}/",
                      forward_files=[item],
                      backward_files=out_list
                      )
        task_list.append(a_task)
        del cmd_list[-1]
    sub_batch_size = None
    if sub_batch_size == None:
        submission = Submission(work_base=path_destination,
                                machine=mach,
                                resources=resrc,
                                task_list=task_list,
                                forward_common_files=[],
                                backward_common_files=[]
                                )
        try:
            submission.run_submission()
        except Exception as e:
            with open('opt_error.log', 'a') as f:
                f.write(str(e) + '\n' + str(task_list) + '\n')

    else:
        num_groups = math.ceil(len(task_list) / sub_batch_size)
        for i in range(num_groups):
            cursor = i * sub_batch_size
            a_task_list = task_list[cursor:cursor + sub_batch_size]
            submission = Submission(work_base=path_destination,
                                    machine=mach,
                                    resources=resrc,
                                    task_list=a_task_list,
                                    forward_common_files=[],
                                    backward_common_files=[]
                                    )
            try:
                submission.run_submission()
            except Exception as e:
                with open('opt_error.log', 'a') as f:
                    f.write(str(e) + '\n' + str(task_list) + '\n')



a=1
max_step = 300
commandlist = [r'/home/mkliu/anaconda3/envs/env001/bin/xtb', '--opt tight']
# commandlist = [r'/home/mkliu/anaconda3/envs/env001/bin/xtb', '--opt tight', '--uhf 1']
out_list = ['xtbopt.xyz', 'xtbopt.log']


run_a_batch(path_source=xyz_path,
            path_destination=dst,
            cmd_list=commandlist,
            out_list=out_list)


os.chdir(dump)
for a_folder in os.listdir(dst):
    if not os.path.isdir(os.path.join(dst,a_folder)):
        continue
    # print(a_folder)
    shutil.copy(src=os.path.join(dst, a_folder, 'xtbopt.xyz'),
                dst=f'{a_folder}.xyz')