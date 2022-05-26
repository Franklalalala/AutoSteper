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
from ase.io.gaussian import read_gaussian_out


hartree2kjmol = Hartree / (kJ / mol)

class Path_Parser():
    def __init__(self, path_para: dict, step: int, start: int, workbase: str, q_cage: Cage, optimizer: Optimizer):
        self.q_cage = q_cage
        self.optimizer = optimizer
        self.workbase = os.path.join(workbase, 'path_info')
        os.makedirs(self.workbase, exist_ok=True)
        self.step = step
        self.start = start
        self.q_path_num = path_para['q_path_num']
        self.q_add_num = path_para['q_add_num']
        self.q_low_e_num = path_para['q_low_e_num']
        self.log_low_e_num = path_para['log_low_e_num']
        if 'ctl_path_para' in path_para.keys():
            self.ctr_path = True
            self.ctl_parent_num = path_para['ctl_path_para']['ctl_parent_num']
            self.max_path_num = path_para['ctl_path_para']['max_path_num']
        else:
            self.ctr_path = False

    def _sort_parent(self, flat_info: pd.DataFrame, info_path: str, parent_flat_info: pd.DataFrame):
        for a_cage in flat_info.keys():
            parent_names = flat_info[a_cage][0]
            name_e_map = {}
            for a_parent in parent_names:
                name_e_map.update({parent_flat_info[a_parent][1]: a_parent})
            new_parent_names = []
            for a_e in sorted(name_e_map.keys())[:self.ctl_parent_num]:
                new_parent_names.append(name_e_map[a_e])
            flat_info[a_cage][0] = new_parent_names
        flat_info.to_pickle(info_path)
        return flat_info



    def get_path_info(self):
        os.chdir(self.workbase)
        def _get_e_area_triangle(e_list: list):
            e_array = np.array(e_list)
            e_array = e_array - self.base_e
            e_area = e_array[1] / 2
            for i in range(2, len(e_array)):
                e_area += (e_array[i] + e_array[i - 1]) / 2
            return e_area * hartree2kjmol

        def _get_e_area(e_list: list):
            e_area = sum(e_list)
            return e_area

        def _get_pathway_unit(idx: int, pathway: list, cage_name: str, e_list: list, name_list: list):
            if idx == 0:
                _, addon_set, _0 = name2seq(cage_size=cage_size, name=cage_name)
                final_pathway = [init_pathway, addon_set, *pathway]
                q_isomer_e = flat_yes_list[idx][cage_name][0]
                final_e_list = [self.base_e, q_isomer_e, *e_list]
                final_name_list = [self.q_cage.name, cage_name, *name_list]

                e_areas.append(_get_e_area(final_e_list))
                pathways.append(final_pathway)
                e_lists.append(final_e_list)
                name_lists.append(final_name_list)
                return
            else:
                new_name_list = [cage_name, *name_list]
                _, addon_set, _0 = name2seq(cage_size=cage_size, name=cage_name)
                parent_name_list = flat_yes_list[idx][cage_name][0]
                q_isomer_e = flat_yes_list[idx][cage_name][1]
                new_e_list = [q_isomer_e, *e_list]
                for a_parent in parent_name_list:
                    _, parent_addon_set, _0 = name2seq(cage_size=cage_size, name=a_parent)
                    new_addon_set = addon_set - parent_addon_set
                    new_pathway = [new_addon_set, *pathway]
                    _get_pathway_unit(idx - 1, pathway=new_pathway, cage_name=a_parent, e_list=new_e_list,
                                      name_list=new_name_list)

        # Get base_e
        if self.optimizer.mode == 'ase':
            atoms = self.q_cage.atoms
            atoms.calc = self.optimizer.calc
            self.base_e = atoms.get_potential_energy()[0][0]
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
        elif self.optimizer.mode == 'gaussian':
            log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
            if not os.path.exists(log_path):
                os.makedirs(name='temp_raw', exist_ok=True)
                os.symlink(src=self.q_cage.pristine_path, dst=os.path.join('temp_raw', self.q_cage.name + '.xyz'))
                os.makedirs(name='temp_opt', exist_ok=True)
                self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt', cmd_list=self.optimizer.cmd_list, cycles=self.optimizer.init_cycle)
            log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
            self.last_image = read_gaussian_out(fd=log_path)[-1]
            self.base_e = self.last_image.get_potential_energy()

        # Dump low e isomers to a log
        deep_info_path = os.path.join(self.q_cage.workbase, f'{self.q_add_num}addons', 'deep_yes_info.pickle')
        deep_info = pd.read_pickle(deep_info_path)
        Max_q_rank = len(deep_info['energy'])
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

        add_num_list = list(range(self.start, self.q_add_num + self.step, self.step))
        flat_yes_list = []

        if self.ctr_path:
            assert len(add_num_list) > 2, print(
                'Current quary addon number is less than 3 steps away from the begining.\n'
                'Do not need control pathes.')
            for ii, i in enumerate(add_num_list):
                info_path = os.path.join(self.q_cage.workbase, f'{i}addons', 'flat_yes_info.pickle')
                a_flat_yes_info = pd.read_pickle(info_path)
                if ii < 2:
                    flat_yes_list.append(a_flat_yes_info)
                else:
                    parent_flat_info = pd.read_pickle(os.path.join(self.q_cage.workbase, f'{i-self.step}addons', 'flat_yes_info.pickle'))
                    sorted_flat_yes_info = self._sort_parent(flat_info=a_flat_yes_info, info_path=info_path, parent_flat_info=parent_flat_info)
                    flat_yes_list.append(sorted_flat_yes_info)
        else:
            for i in add_num_list:
                info_path = os.path.join(self.q_cage.workbase, f'{i}addons', 'flat_yes_info.pickle')
                a_flat_yes_info = pd.read_pickle(info_path)
                flat_yes_list.append(a_flat_yes_info)

        # Query the path of some low e isomers
        self.q_low_e_num = min(self.q_low_e_num, Max_q_rank)
        for q_rank in range(self.q_low_e_num):
            q_isomer_e = deep_info['energy'][q_rank]
            q_isomer_name = deep_info['name'][q_rank]
            q_isomer_atoms = read(deep_info['xyz_path'][q_rank])
            q_rank_workbase = os.path.join(self.workbase, f'isomer_rank_{q_rank}')
            os.makedirs(exist_ok=True, name=q_rank_workbase)
            os.chdir(q_rank_workbase)
            write(filename=q_isomer_name + '.xyz', images=q_isomer_atoms, format='xyz', comment=str(q_isomer_e))
            init_pathway = {0}

            idx = len(add_num_list) - 1
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
            if self.ctr_path:
                sorted_info = sorted_info[sorted_info.index < self.max_path_num]
            sorted_info.to_pickle(info_path)

            # Get relatively energy info
            e_lists = sorted_info['e_list']
            e_array = np.array(e_lists[0])
            for i in e_lists[1:]:
                e_array = np.vstack((e_array, np.array(i)))

            if len(e_array.shape) == 1:
                print(f'There are only one root for rank {q_rank}.')
                continue

            rel_e_array = np.ones_like(e_array)
            for i in range(len(add_num_list)+1):
                rel_e_array[:, i] = e_array[:, i] - min(e_array[:, i])
            rel_e_array = rel_e_array * hartree2kjmol

            # Dump path related xyz file
            for path_rank in range(min(self.q_path_num, len(e_lists))):
                name_list = sorted_info['name'][path_rank]
                one_path_workbase = os.path.join(q_rank_workbase, f'path_rank_{path_rank}')
                os.makedirs(exist_ok=True, name=one_path_workbase)
                os.chdir(one_path_workbase)
                log_path = f'path_rank_{path_rank}_root.log'
                write(filename=log_path, images=self.q_cage.atoms, format='xyz', append=True)
                for idx, add_num in enumerate(add_num_list):
                    addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                    name = name_list[idx+1]
                    new_path = f'{name}_addon_{add_num}.xyz'
                    if self.optimizer.mode == 'xtb':
                        xyz_filename = 'xtbopt.xyz'
                    elif self.optimizer.mode == 'ase':
                        xyz_filename = 'opt.xyz'
                    elif self.optimizer.mode == 'gaussian':
                        xyz_filename = 'gau.xyz'

                    opt_path = os.path.join(addon_path, f'opt_final', name, xyz_filename)
                    if os.path.exists(opt_path):
                        shutil.copy(src=opt_path, dst=new_path)
                    else:
                        opt_path = os.path.join(addon_path, f'opt_{self.optimizer.init_cycle}', name, xyz_filename)
                        shutil.copy(src=opt_path, dst=new_path)
                    atoms = read(new_path, format='xyz')
                    write(filename=log_path, images=atoms, format='xyz', append=True)

            np.save(file=os.path.join(q_rank_workbase, f'Path_relative_energy.npy'), arr=rel_e_array)
            # Plot relative energy for path rank
            rel_e_array = rel_e_array[:self.q_path_num, :]
            fig_1 = plt.figure(dpi=400, figsize=(len(add_num_list), self.q_path_num))
            cmap = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
            rel_e_df = pd.DataFrame(rel_e_array)
            rel_e_df.columns = [0, *add_num_list]
            sns.heatmap(rel_e_df, annot=True, cmap=cmap, linewidths=.5)
            plt.ylabel('Path rank.')
            plt.xlabel('Addon number.')
            plt.title('Path relative energy (kj/mol).')
            plt.savefig(os.path.join(q_rank_workbase, f'Path_relative_energy.png'))

            # Plot relative energy for isomers.
            isomer_rel_e = np.zeros_like(e_array)
            for idx, add_num in enumerate(add_num_list):
                deep_info_path = os.path.join(self.q_cage.workbase, f'{add_num}addons', 'deep_yes_info.pickle')
                deep_info = pd.read_pickle(deep_info_path)
                isomer_min = deep_info['energy'][0]
                isomer_rel_e[:, idx+1] = e_array[:, idx+1] - isomer_min
            isomer_rel_e = isomer_rel_e * hartree2kjmol


            np.save(file=os.path.join(q_rank_workbase, f'Isomer_relative_energy.npy'), arr=isomer_rel_e)
            isomer_rel_e = isomer_rel_e[:self.q_path_num, :]
            # Simple plot
            fig_2 = plt.figure(dpi=400, figsize=(len(add_num_list), self.q_path_num))
            cmap = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
            isomer_rel_e_df = pd.DataFrame(isomer_rel_e)
            isomer_rel_e_df.columns = [0, *add_num_list]
            sns.heatmap(isomer_rel_e_df, annot=True, cmap=cmap, linewidths=.5)
            plt.ylabel('Path rank.')
            plt.xlabel('Addon number.')
            plt.title('Isomer relative energy (kj/mol).')
            plt.savefig(os.path.join(q_rank_workbase, f'Isomer_relative_energy.png'))
