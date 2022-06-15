import os
import shutil

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import tqdm
from ase import Atoms
from ase.io import read, write
from ase.io.gaussian import read_gaussian_out
from ase.neighborlist import build_neighbor_list
from ase.units import Hartree, kJ, mol
from autosteper.cage import name2seq, Cage, seq2name
from autosteper.optimizers import *
from autosteper.tools import strip_extraFullerene
from tqdm import tqdm


def get_G(atoms: Atoms):
    nb_list = build_neighbor_list(atoms)
    mat = nb_list.get_connectivity_matrix()
    adj = mat.toarray()
    G = nx.from_numpy_array(adj)
    return G


def test_iso(q_atoms: Atoms, deep_info_path: str, top_k: int=None):
    print('Please make sure the queried atoms has been optimized with at least semi-empirical level of theory.')
    target_G = get_G(q_atoms)
    deep_info = pd.read_pickle(deep_info_path)
    flag = 0
    for idx, a_xyz_path in enumerate(tqdm(deep_info['xyz_path'][:top_k])):
        an_atoms = read(a_xyz_path)
        pred_G = get_G(an_atoms)
        if nx.is_isomorphic(pred_G, target_G):
            rel_e = (deep_info['energy'][idx] - deep_info['energy'][0])*Hartree
            print('The queried atoms are found in AutoStepper scan results.\n'
                  f'The target is in rank {str(idx+1)}.\n'
                  f'Has a relative energy of {str(rel_e)}eV considering the lowest isomer in scan results.\n')
            flag = 1
            break
    if flag == 0:
        if top_k:
            print(f'The queried atoms are NOT found in top {top_k} AutoStepper scan results.\n'
                  'Please check its substructure with test_sub function or enlarge test number.')
        else:
            print('The queried atoms are NOT found in AutoStepper scan results.\n'
                  'Please check its substructure with test_sub function.')
        return None
    else:
        return idx+1


def test_sub(workbase: str, dump_file_path: str, step: int, first_addon: int, cage_size: int=None,q_atoms: Atoms=None, q_seq: str=None, q_cage: Cage=None, q_add_set: set=None):
    if q_atoms:
        print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
              'The atom sequence in two cages should be identical.\n'
              'The addons should be atomic addons.\n'
              'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
        _, q_add_set = strip_extraFullerene(q_atoms)
    elif q_seq:
        _, q_add_set, _0 = seq2name(q_seq, q_cage)
    elif not q_add_set:
        print('Please input query addon set.\n'
              'Currently only surpport ase.atoms input, sequence input and addon set input.')
    max_add_num = len(q_add_set)
    if cage_size == None:
        if q_cage:
            cage_size = q_cage.size
        else:
            print('Please input cage size.')

    for an_addon in range(first_addon, max_add_num+step, step):
        addon_path = str(an_addon) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)
        with open(dump_file_path, 'a') as f:
            deep_info = pd.read_pickle('deep_yes_info.pickle')
            min_e = deep_info['energy'][0]
            f.write(f'Start screen substructures in {addon_path}:\n')
            flag = 0
            for idx, a_name in enumerate(deep_info['name']):
                _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
                if len(pred_add_set - q_add_set) == 0:
                    f.write('A substructrue is found:\n')
                    f.write(f'The absent addon sites are:      {str(q_add_set - pred_add_set)}.\n')
                    f.write(f'In rank:                         {str(idx+1)}\n')
                    rel_e = (deep_info['energy'][idx] - min_e)*Hartree
                    f.write(f'Relative energy is:              {str(rel_e)}eV\n')
                    f.write(f'Name is:                         {a_name}\n')
                    f.write('\n')
                    f.write('\n')
                    flag = 1
                    # break
            if flag == 0:
                f.write(f'{addon_path} do not have a substructure.\n')

def test_sub_plot(workbase: str, dump_pic_path: str, step: int, first_addon: int, cage_size: int=None,q_atoms: Atoms=None, q_seq: str=None, q_cage: Cage=None, q_add_set: set=None, cut_rank: int=None, cut_value: float=None):
    if q_atoms:
        print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
              'The atom sequence in two cages should be identical.\n'
              'The addons should be atomic addons.\n'
              'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
        _, q_add_set = strip_extraFullerene(q_atoms)
    elif q_seq:
        _, q_add_set, _0 = seq2name(q_seq, q_cage)
    elif not q_add_set:
        print('Please input query addon set.\n'
              'Currently only surpport ase.atoms input, sequence input and addon set input.')
    max_add_num = len(q_add_set)
    if cage_size == None:
        if q_cage:
            cage_size = q_cage.size
        else:
            print('Please input cage size.')

    add_list = list(range(first_addon, max_add_num+step, step))
    fig = plt.figure(dpi=400, figsize=(len(add_list), 7))
    is_first = True
    seed_info_dict = dict()
    non_seed_info_dict = dict()
    for an_addon in add_list:
        addon_path = str(an_addon) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)

        deep_info = pd.read_pickle('deep_yes_info.pickle')
        min_e = deep_info['energy'][0]

        rel_seed = []
        rel_non_seed = []
        old_rel_seed = []
        old_rel_non_seed = []

        for idx, a_name in enumerate(deep_info['name']):
            _, pred_add_set, _0 = name2seq(name=a_name, cage_size=84)
            if len(pred_add_set - q_add_set) == 0:
                rel_e = (deep_info['energy'][idx] - min_e) * Hartree
                old_rel_e = (deep_info['energy'][idx]) * Hartree
                if idx < cut_rank and rel_e < cut_value:
                    old_rel_seed.append(old_rel_e)
                    rel_seed.append(old_rel_e)
                else:
                    old_rel_non_seed.append(old_rel_e)
                    rel_non_seed.append(old_rel_e)

        if len(rel_seed)+len(rel_non_seed)>0:
            if is_first:
                plt.scatter(x=[an_addon] * len(rel_seed), y=rel_seed, marker='+', c='red', label='seed')
                plt.scatter(x=[an_addon] * len(rel_non_seed), y=rel_non_seed, marker='x', c='blue', label='non-seed')
                is_first = False
            else:
                plt.scatter(x=[an_addon] * len(rel_seed), y=rel_seed, marker='+', c='red')
                plt.scatter(x=[an_addon] * len(rel_non_seed), y=rel_non_seed, marker='x', c='blue')
        seed_info_dict.update({an_addon: [old_rel_seed]})
        non_seed_info_dict.update({an_addon: [old_rel_non_seed]})
    seed_info_df = pd.DataFrame(seed_info_dict)
    non_seed_info_df = pd.DataFrame(non_seed_info_dict)
    seed_info_df.to_pickle(dump_pic_path[:-4]+'_seed.pickle')
    non_seed_info_df.to_pickle(dump_pic_path[:-4]+'_nonseed.pickle')

    plt.legend(loc='best')
    plt.ylim(None,1)
    plt.xticks(ticks=add_list, labels=add_list)
    plt.xlabel('Addon number')
    plt.ylabel('Substructure relative energy (ev).')
    # plt.title(f'{dump_pic_path[:-4]}')
    plt.savefig(dump_pic_path)


hartree2kjmol = Hartree / (kJ / mol)
class Path_Parser():
    def __init__(self, path_para: dict, step: int, start: int, workbase: str, q_cage: Cage, optimizer: Optimizer, refine_para: dict=None):
        self.q_cage = q_cage
        self.optimizer = optimizer
        self.workbase = os.path.join(workbase, f'path_info_{self.q_cage.name}')
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

        if 'refine_para' in path_para.keys():
            self.is_refine = True
            self.refine_top_k = path_para['refine_para']['refine_top_k']
            self.refine_mode = path_para['refine_para']['opt_mode']
            if self.refine_mode == 'xtb':
                self.refine_optimizer = XTB_Optimizer(path_para['refine_para']['opt_para'], checker=None, cage=None)
            elif self.refine_mode == 'ase':
                self.refine_optimizer = ASE_Optimizer(path_para['refine_para']['opt_para'], checker=None, cage=None)
            elif self.refine_mode == 'gaussian':
                self.refine_optimizer = Gaussian_Optimizer(path_para['refine_para']['opt_para'], checker=None, cage=None)
        else:
            self.is_refine = False


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

        def _get_base_e():
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
                    self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt',
                                               cmd_list=self.optimizer.cmd_list)
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
                    self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt',
                                               cmd_list=self.optimizer.cmd_list, cycles=self.optimizer.init_cycle)
                log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
                self.last_image = read_gaussian_out(fd=log_path)[-1]
                self.base_e = self.last_image.get_potential_energy()

        def _refine_e(old_info: pd.DataFrame):
            if self.refine_optimizer.mode == 'ase':
                atoms = self.q_cage.atoms
                atoms.calc = self.optimizer.calc
                self.base_e = atoms.get_potential_energy()[0][0]
                return None
            elif self.refine_optimizer.mode == 'xtb':
                os.makedirs(name='refine_workbase', exist_ok=True)
                os.chdir('refine_workbase')
                os.makedirs(name='odd_raw', exist_ok=True)
                os.makedirs(name='even_raw', exist_ok=True)
                os.makedirs(name='opt_workbase', exist_ok=True)

                for add_num_idx, add_num in enumerate(add_num_list[:-1]):
                    addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                    for name_list in old_info['name']:
                        name = name_list[add_num_idx + 1]
                        if add_num % 2 == 1:
                            new_path = os.path.join('odd_raw', f'{name}.xyz')
                        else:
                            new_path = os.path.join('even_raw', f'{name}.xyz')
                        if os.path.exists(new_path):
                            continue
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
                self.refine_optimizer.run_a_batch(path_source='even_raw', path_destination='opt_workbase',
                                                  cmd_list=self.refine_optimizer.cmd_list)
                self.refine_optimizer.run_a_batch(path_source='odd_raw', path_destination='opt_workbase',
                                                  cmd_list=[*self.refine_optimizer.cmd_list, '--uhf 1'])
                os.chdir('opt_workbase')
                name_e_map = dict()
                for a_name in os.listdir('./'):
                    with open(os.path.join(a_name, 'xtbopt.xyz'), 'r') as f:
                        f.readline()
                        energy_line = f.readline()
                        energy = float(energy_line.split()[1])
                    name_e_map.update({a_name: energy})

                for add_num_idx, add_num in enumerate(add_num_list[:-1]):
                    addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                    for pathway_idx, name_list in enumerate(old_info['name']):
                        name = name_list[add_num_idx + 1]
                        e = name_e_map[name]
                        old_info['e_list'][pathway_idx][add_num_idx + 1] = e
                for idx, e_list in enumerate(old_info['e_list']):
                    e_area = _get_e_area(e_list=e_list)
                    old_info['e_area'][idx] = e_area

                os.chdir(self.workbase)
                refined_info = old_info.sort_values(by='e_area')
                refined_info.index = sorted(refined_info.index)
                return refined_info
            # elif self.optimizer.mode == 'gaussian':
            #     log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
            #     if not os.path.exists(log_path):
            #         os.makedirs(name='temp_raw', exist_ok=True)
            #         os.symlink(src=self.q_cage.pristine_path, dst=os.path.join('temp_raw', self.q_cage.name + '.xyz'))
            #         os.makedirs(name='temp_opt', exist_ok=True)
            #         self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt',
            #                                    cmd_list=self.optimizer.cmd_list, cycles=self.optimizer.init_cycle)
            #     log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
            #     self.last_image = read_gaussian_out(fd=log_path)[-1]
            #     self.base_e = self.last_image.get_potential_energy()

        self.base_e = 0

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
            info = pd.DataFrame({'pathway': pathways, 'name': name_lists, 'e_list': e_lists, 'e_area': e_areas})
            sorted_info = info.sort_values(by='e_area')
            sorted_info.index = sorted(sorted_info.index)

            if self.ctr_path:
                sorted_info = sorted_info[sorted_info.index < self.max_path_num]
            if self.is_refine:
                sorted_info = _refine_e(old_info=sorted_info[sorted_info.index < self.refine_top_k])

            info_path = f'{q_isomer_name}_pathway_info.pickle'
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
            fig_1 = plt.figure(dpi=400, figsize=(len(add_num_list), min(self.q_path_num, len(e_lists) + 1, 100)))
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
            fig_2 = plt.figure(dpi=400, figsize=(len(add_num_list), min(self.q_path_num, len(e_lists) + 1, 100)))
            cmap = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
            isomer_rel_e_df = pd.DataFrame(isomer_rel_e)
            isomer_rel_e_df.columns = [0, *add_num_list]
            sns.heatmap(isomer_rel_e_df, annot=True, cmap=cmap, linewidths=.5)
            plt.ylabel('Path rank.')
            plt.xlabel('Addon number.')
            plt.title('Isomer relative energy (kj/mol).')
            plt.savefig(os.path.join(q_rank_workbase, f'Isomer_relative_energy.png'))
