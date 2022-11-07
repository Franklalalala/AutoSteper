import os
import shutil

import numpy as np
import pandas as pd
from ase.atoms import Atoms
from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list
from autosteper.cage import name2seq, Cage


class Checker():
    def __init__(self, pst_cage: Cage, group: str, skin: float = 0.15):
        self.skin = 0.15
        self.group = group
        self.pst_cage = pst_cage

    def check_ADJ(self, dok_matrix):
        for idx, an_atom in enumerate(self.q_atoms):
            dok_adj = dok_matrix[idx]
            adj_array = dok_adj.tocoo().col
            len_adj_arr = len(adj_array)
            if an_atom.symbol == self.pst_cage.symbol:
                if len_adj_arr == 4:
                    flag = 0
                    for ii in adj_array:
                        if self.q_atoms[ii].symbol in [self.group, 'H', 'O']:
                            flag = 1
                            break
                    if flag:
                        continue
                    else:
                        return 3
                elif len_adj_arr < 3:
                    return 4
                elif len_adj_arr == 3:
                    for ii in adj_array:
                        if self.q_atoms[ii].symbol in [self.group, 'H', 'O']:
                            return 4
                elif len_adj_arr > 4:
                    return 6
            elif len_adj_arr >= 2:
                return 5
            elif len_adj_arr == 0:
                return 1
        return 0

    def check_group_intactness(self):
        dummy_image = Atoms()
        neighborlist = build_neighbor_list(self.q_atoms, bothways=True, self_interaction=False, skin=self.skin)
        dok_matrix = neighborlist.get_connectivity_matrix()
        if self.group == 'OH':
            for idx, an_atom in enumerate(self.q_atoms):
                if an_atom.symbol == 'H':
                    dok_adj = dok_matrix[idx]
                    adj_array = dok_adj.tocoo().col
                    if not len(adj_array) == 1:
                        return 7, None
                else:
                    dummy_image.append(an_atom)

        if self.group == 'CH3':
            for idx, an_atom in enumerate(self.q_atoms):
                if an_atom.symbol == 'H':
                    dok_adj = dok_matrix[idx]
                    adj_array = dok_adj.tocoo().col
                    if not len(adj_array) == 1:
                        return 7, None
                elif idx > self.pst_cage.size - 1:
                    an_atom.symbol = 'H'
                    dummy_image.append(an_atom)
                else:
                    dummy_image.append(an_atom)

        if self.group == 'CF3':
            for idx, an_atom in enumerate(self.q_atoms):
                if an_atom.symbol == 'F':
                    dok_adj = dok_matrix[idx]
                    adj_array = dok_adj.tocoo().col
                    if not len(adj_array) == 1:
                        return 7, None
                elif idx > self.pst_cage.size - 1:
                    an_atom.symbol = 'H'
                    dummy_image.append(an_atom)
                else:
                    dummy_image.append(an_atom)
        return None, dummy_image

    def check_a_job(self, q_atoms: Atoms, q_name: str):
        _, init_addonset, bin_arr = name2seq(q_name, cage_size=self.pst_cage.size)
        self.q_atoms = q_atoms
        if self.group in ['CH3', 'CF3', 'OH']:
            group_status, self.q_atoms = self.check_group_intactness()
            if group_status:
                return group_status
        if self.group in ['OH', 'Cl', 'F', 'H', 'Br']:
            cutoffs = natural_cutoffs(self.q_atoms)
        else:
            cutoffs = [natural_cutoffs(self.pst_cage.atoms)[0]] * len(self.q_atoms)
        neighborlist = NeighborList(cutoffs, skin=self.skin, bothways=True, self_interaction=False)
        neighborlist.update(self.q_atoms)
        dok_matrix = neighborlist.get_connectivity_matrix()
        ADJ_status = self.check_ADJ(dok_matrix=dok_matrix)
        if ADJ_status:
            return ADJ_status
        bonds = np.array([*dok_matrix.keys()])
        map_dict = dict()
        for i in range(bonds.shape[0]):
            bond = bonds[i]
            symbol_0 = self.q_atoms[bond[0]].symbol
            symbol_1 = self.q_atoms[bond[1]].symbol
            if not symbol_0 == symbol_1:
                if bond[0] < bond[1]:
                    if symbol_0 == self.pst_cage.symbol:
                        cage_element_idx = bond[0]
                        addon_idx = bond[1]
                    else:
                        cage_element_idx = bond[1]
                        addon_idx = bond[0]
                    map_dict.update({addon_idx: cage_element_idx})
        opted_addonset = set(map_dict.values())
        if init_addonset == opted_addonset:
            return 0
        elif len(init_addonset) == len(opted_addonset):
            return 2
        else:
            return 1

    # job_status == 0 means the topology stays unchanged during optimization.
    # job_status == 1 means at least one addon atom breaks the bond with the cage and becomes a radical during optimization.
    # job_status == 2 means at least one addon atom breaks the bond with the original addon site and changed to another during optimization.
    # job_status == 3 means at least one 3 membered carbon ring formed during optimization, which is extremely unstable. (carbon for instance, this can be changed for other elements.)
    # job_status == 4 means at least one carbon atom only has 2 neighbor carbons or less, which means the cage is broken.
    # job_status == 5 means at least one addon atom binds with 2 or more atoms.
    # job_status == 6 means at least one carbon atom has 5 or more neighbors, which means a small cluster is formed.
    # job_status == 7 means the inner intactness of at least one addon group has broken.
    def check(self, passed_info_path: str, status_info_path: str):
        check_cwd_ = os.getcwd()
        passed_info = pd.read_pickle(passed_info_path)
        status_info = pd.read_pickle(status_info_path)
        failed_idxes = []
        failed_status_codes = []
        failed_job_paths = []
        q_name = os.path.splitext(os.path.basename(passed_info['xyz_path'][0]))[0]
        failed_folder_path = os.path.join(check_cwd_, 'cooking')
        if q_name in os.listdir(failed_folder_path):
            for idx, a_xyz_path in enumerate(passed_info['xyz_path']):
                q_atoms = read(a_xyz_path)
                q_name = os.path.splitext(os.path.basename(a_xyz_path))[0]
                job_status = self.check_a_job(q_atoms=q_atoms, q_name=q_name)
                if job_status == 0:
                    continue
                status_info[q_name] = job_status
                failed_status_codes.append(job_status)
                failed_idxes.append(idx)
                failed_job_paths.append(os.path.join(failed_folder_path, q_name))
                os.remove(path=a_xyz_path)
        else:
            failed_folder_path = os.path.join(check_cwd_, 'failed')
            os.makedirs(failed_folder_path, exist_ok=True)
            for idx, a_xyz_path in enumerate(passed_info['xyz_path']):
                q_atoms = read(a_xyz_path)
                q_name = os.path.splitext(os.path.basename(a_xyz_path))[0]
                job_status = self.check_a_job(q_atoms=q_atoms, q_name=q_name)
                if job_status == 0:
                    continue
                status_info[q_name] = job_status
                failed_status_codes.append(job_status)
                failed_idxes.append(idx)
                failed_new_path = os.path.join(failed_folder_path, q_name + '.xyz')
                shutil.move(src=a_xyz_path, dst=failed_new_path)
                failed_job_paths.append(failed_new_path)

        passed_info = passed_info.drop(failed_idxes)
        passed_info.index = list(range(len(passed_info.index)))
        passed_info.to_pickle(passed_info_path)
        status_info.to_pickle(status_info_path)
        status_codes = list(status_info.T[0])

        failed_info = ''
        for idx, failed_idx in enumerate(failed_status_codes):
            failed_info = failed_info + f'{failed_idx}        {failed_job_paths[idx]}\n'
        with open('failed_job_paths', 'w') as f:
            f.write(failed_info)
        return status_codes
