import os
import shutil
from pathlib import Path
from typing import Union
import networkx as nx
from ase.neighborlist import build_neighbor_list
from ase.visualize import view
from tqdm import tqdm
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
from autosteper.optimizers import Optimizer
from autosteper.cage import name2seq, Cage
from ase.units import Hartree, kJ, mol
import numpy as np
import pandas as pd
from ase import Atoms, Atom
from ase.io import read, write


hartree2kjmol = Hartree / (kJ / mol)

class Path_Parser():
    def __init__(self, path_para: dict, step: int, workbase: str, q_cage: Cage, optimizer: Optimizer):
        self.q_cage = q_cage
        self.optimizer = optimizer
        self.workbase = os.path.join(workbase, 'path_info')
        os.makedirs(self.workbase, exist_ok=True)
        self.step = step
        self.q_path_num = path_para['q_path_num']
        self.q_add_num = path_para['q_add_num']
        self.q_low_e_num = path_para['q_low_e_num']
        self.log_low_e_num = path_para['log_low_e_num']



    def get_path_info(self):
        os.chdir(self.workbase)
        def _get_e_area(e_list: list):
            e_array = np.array(e_list)
            for i in range(len(e_array)):
                e_array = e_array - self.base_e
            e_area = e_array[0] / 2
            for i in range(1, len(e_array)):
                e_area += (e_array[i] + e_array[i - 1]) / 2
            return e_area

        def _get_pathway_unit(idx: int, pathway: list, cage_name: str, e_list: list, name_list: list):
            if idx == 0:
                _, addon_set = name2seq(cage_size=cage_size, name=cage_name)
                final_pathway = [init_pathway, addon_set, *pathway]
                q_isomer_e = flat_yes_list[idx][cage_name][0]
                final_e_list = [self.base_e, q_isomer_e, *e_list]
                final_name_list = [self.q_cage.pristine_path, cage_name, *name_list]

                e_areas.append(_get_e_area(final_e_list))
                pathways.append(final_pathway)
                e_lists.append(final_e_list)
                name_lists.append(final_name_list)
                return
            else:
                new_name_list = [cage_name, *name_list]
                _, addon_set = name2seq(cage_size=cage_size, name=cage_name)
                parent_name_list = flat_yes_list[idx][cage_name][0]
                q_isomer_e = flat_yes_list[idx][cage_name][1]
                new_e_list = [q_isomer_e, *e_list]
                for a_parent in parent_name_list:
                    _, parent_addon_set = name2seq(cage_size=cage_size, name=a_parent)
                    new_addon_set = addon_set - parent_addon_set
                    new_pathway = [new_addon_set, *pathway]
                    _get_pathway_unit(idx - 1, pathway=new_pathway, cage_name=a_parent, e_list=new_e_list,
                                      name_list=new_name_list)

        # Get base_e
        if self.optimizer.mode == 'ase':
            pass
        elif self.optimizer.mode == 'xtb':
            opted_xyz = os.path.join('temp_opt', self.q_cage.name, 'xtbopt.xyz')
            if not os.path.exists(opted_xyz):
                os.makedirs(name='temp_raw', exist_ok=True)
                os.symlink(src=self.q_cage.pristine_path, dst=os.path.join('temp_raw', self.q_cage.name + '.xyz'))
                os.makedirs(name='temp_opt', exist_ok=True)
                self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt', cmd_list=self.optimizer.cmd_list)
            with open(opted_xyz, 'r') as f:
                f.readline()
                e_line = f.readline()
                self.base_e = float(e_line.split()[1])

        # Dump low e isomers to a log
        deep_info_path = os.path.join(self.q_cage.workbase, f'{self.q_add_num}addons', 'deep_yes_info.pickle')
        deep_info = pd.read_pickle(deep_info_path)
        Max_q_rank = len(deep_info['energy'])
        self.q_low_e_num = min(self.q_low_e_num, Max_q_rank)
        self.log_low_e_num = min(self.log_low_e_num, Max_q_rank)

        log_path = os.path.join(f'top_{self.log_low_e_num}_isomers.log')
        for rank in range(self.log_low_e_num):
            src_path = deep_info['xyz_path'][rank]
            if rank == 0:
                shutil.copy(src=src_path, dst=log_path)
            else:
                with open(log_path, 'a') as w_file, open(src_path, 'r') as r_file:
                    for a_line in r_file.readlines():
                        w_file.write(a_line)


        # Query the path of some low e isomers
        for q_rank in range(self.q_low_e_num):
            q_isomer_e = deep_info['energy'][q_rank]
            q_isomer_name = deep_info['name'][q_rank]
            q_isomer_atoms = read(deep_info['xyz_path'][q_rank])
            q_rank_workbase = os.path.join(self.workbase, f'rank_{q_rank}')
            os.makedirs(exist_ok=True, name=q_rank_workbase)
            os.chdir(q_rank_workbase)
            write(filename=q_isomer_name + '.xyz', images=q_isomer_atoms, format='xyz', comment=str(q_isomer_e))
            init_pathway = {0}
            flat_yes_list = []
            for i in range(self.step, self.q_add_num + self.step, self.step):
                addon_path = os.path.join(self.q_cage.workbase, f'{i}addons')
                a_flat_yes_info = pd.read_pickle(os.path.join(addon_path, 'flat_yes_info.pickle'))
                flat_yes_list.append(a_flat_yes_info)

            idx = int(self.q_add_num / self.step) - 1
            pathways = []
            e_lists = []
            e_areas = []
            name_lists = []
            cage_size = self.q_cage.size

            # Get path info recursively
            _get_pathway_unit(idx=idx, pathway=[], cage_name=q_isomer_name, e_list=[], name_list=[])
            info_path = f'{q_isomer_name}_pathway_info.pickle'
            info = pd.DataFrame({'pathway': pathways, 'name': name_lists, 'e_list': e_lists, 'e_area': e_areas})
            sorted_info = info.sort_values(by='e_area')
            sorted_info.index = sorted(sorted_info.index)
            sorted_info.to_pickle(info_path)

            # Get relatively energy info
            e_array = np.array(e_lists[0])
            for i in e_lists[1:]:
                e_array = np.vstack((e_array, np.array(i)))
            rel_e_array = np.ones_like(e_array)

            if len(e_array.shape) == 1:
                print(f'There are only one root for rank {q_rank}.')
                continue

            num_steps = e_array.shape[-1]
            for i in range(num_steps):
                rel_e_array[:, i] = e_array[:, i] - min(e_array[:, i])

            rel_e_array = rel_e_array * hartree2kjmol

            self.q_path_num = min(self.q_path_num, len(e_lists))
            # Dump path related xyz file
            for path_rank in range(self.q_path_num):
                name_list = sorted_info['name'][path_rank]
                one_path_workbase = os.path.join(q_rank_workbase, str(path_rank))
                os.makedirs(exist_ok=True, name=one_path_workbase)
                os.chdir(one_path_workbase)
                log_path = f'rank_{path_rank}_root.log'
                write(filename=log_path, images=self.q_cage.atoms, format='xyz', append=True)
                for idx, name in enumerate(name_list):
                    if idx == 0:
                        continue
                    addon_path = os.path.join(self.q_cage.workbase, f'{self.step * idx}addons')
                    new_path = f'{name}_step_{self.step * (idx + 1) - 1}.xyz'
                    if self.optimizer.mode == 'xtb':
                        xyz_filename = 'xtbopt.xyz'
                    elif self.optimizer.mode == 'ase':
                        xyz_filename = 'opt.xyz'

                    opt_path = os.path.join(addon_path, f'opt_final', name, xyz_filename)
                    if os.path.exists(opt_path):
                        shutil.copy(src=opt_path, dst=new_path)
                    else:
                        opt_path = os.path.join(addon_path, f'opt_{self.optimizer.init_cycle}', name, xyz_filename)
                        shutil.copy(src=opt_path, dst=new_path)
                    atoms = read(new_path, format='xyz')
                    write(filename=log_path, images=atoms, format='xyz', append=True)

            # Simple plot
            x = np.linspace(start=0, stop=num_steps - 1, num=num_steps)
            fig = plt.figure()
            for i in range(self.q_path_num):
                plt.plot(rel_e_array[i], label=f'rank_{i}')
                plt.scatter(x=x, y=rel_e_array[i], label=f'rank_{i}', marker='+')
            plt.ylabel('Relative energy (kj/mol)')
            plt.xlabel('Addon number.')
            plt.legend(loc='best')
            plt.savefig(os.path.join(q_rank_workbase, f'relative_energy_rank_{q_rank}.png'))
            np.save(file=os.path.join(q_rank_workbase, f'relative_energy_rank_{q_rank}.npy'), arr=rel_e_array)


